/**
 * Author: Kun Sun (sunkun@szbl.ac.cn)
 * Date: Feb, 2020
 * This program is part of the Ktrim package
**/

#include <fstream>
#include <sstream>
#include <algorithm>
#include <thread>
#include <ctime>
#include <stdlib.h>
#include <memory.h>
#include <omp.h>
#include "common.h"
#include <functional>
#include <mutex>
using namespace std;

typedef rabbit::core::TDataQueue<rabbit::fq::FastqDataPairChunk> FqDataPairChunkQueue;
std::mutex buffer_head_mtx;
std::mutex buffer_size_mtx;

void inline CPEREAD_resize( CPEREAD * read, int n ) {
	read->size = n;
	read->seq1[  n ] = 0;
	read->qual1[ n ] = 0;
	read->seq2[  n ] = 0;
	read->qual2[ n ] = 0;
}

bool inline is_revcomp( const char a, const char b ) {
/*
	if( a=='A' )
		return b=='T';
	if( a=='C' )
		return b=='G';
	if( a=='G' )
		return b=='C';
	if( a=='T' )
		return b=='A';

	//TODO: consider how to deal with N, call it positive or negative???
	return false;
*/
	switch( a ) {
		case 'A': return b=='T';
		case 'C': return b=='G';
		case 'G': return b=='C';
		case 'T': return b=='A';
		default : return false;
	}
}

void find_seed_pe( vector<unsigned int> &seed, const CPEREAD *read, const ktrim_param & kp ) {
	seed.clear();
	register const char *poffset  = read->seq1;
	register const char *indexloc = poffset;
	while( true ) {
		indexloc = strstr( indexloc, kp.adapter_index1 );
		if( indexloc == NULL )
			break;
		seed.push_back( indexloc - poffset );
		++ indexloc;
	}
	poffset  = read->seq2;
	indexloc = poffset;
	while( true ) {
		indexloc = strstr( indexloc, kp.adapter_index2 );
		if( indexloc == NULL )
			break;
		seed.push_back( indexloc - poffset );
		++ indexloc;
	}
	sort( seed.begin(), seed.end() );
}

// this function is slower than C++ version
void workingThread_PE_C( unsigned int tn, unsigned int start, unsigned int end, CPEREAD *workingReads,
					ktrim_stat * kstat, writeBufferTotal * write_buffer_total, const ktrim_param & kp ) {

	int b1stored  = 0, b2stored = 0;

//	vector<unsigned int> seed;
//	vector<unsigned int> :: iterator it, end_of_seed;
	register int *seed = new int[ MAX_SEED_NUM ];
	register int hit_seed;
	register int *it, *end_of_seed;

	register CPEREAD *wkr = workingReads;
	for( unsigned int iii=end-start; iii; --iii, ++wkr ) {
		// quality control
		register int i = get_quality_trim_cycle_pe( wkr, kp ); // i 是进行qulity-trim 的位置
		if( i == 0 ) { // not long enough
			++ kstat->dropped[ tn ];
			continue;
		}
		if( i != wkr->size ) {  // quality-trim occurs
			CPEREAD_resize( wkr, i );
		}
 
		hit_seed = 0;
		register const char *poffset  = wkr->seq1;
		register const char *indexloc = poffset;
		while( true ) {
			indexloc = strstr( indexloc, kp.adapter_index1 );
			if( indexloc == NULL )
				break;
			//seed.push_back( indexloc - poffset );
			seed[ hit_seed++ ] = indexloc - poffset;
			++ indexloc;
		}
		poffset  = wkr->seq2;
		indexloc = poffset;
		while( true ) {
			indexloc = strstr( indexloc, kp.adapter_index2 );
			if( indexloc == NULL )
				break;
			//seed.push_back( indexloc - poffset );
			seed[ hit_seed++ ] = indexloc - poffset;
			++ indexloc;
		}
		//sort( seed.begin(), seed.end() );
		end_of_seed = seed + hit_seed;
		if( hit_seed != 0 )
			sort( seed, seed + hit_seed );

		register unsigned int last_seed = impossible_seed;	// a position which cannot be a seed
		//end_of_seed = seed.end();
		//for( it=seed.begin(); it!=end_of_seed; ++it ) {
		for( it=seed; it!=end_of_seed; ++it ) {
			if( *it != last_seed ) {
			// as there maybe the same value in seq1_seed and seq2_seed,
			// use this to avoid re-calculate that pos
				if( check_mismatch_dynamic_PE_C( wkr, *it, kp ) )
					break;

				last_seed = *it;
			}
		}
		if( it != end_of_seed ) {	// adapter found
			++ kstat->real_adapter[tn];
			if( (*it)+1 >= kp.min_length )	{
				CPEREAD_resize( wkr, (*it)+1 );
			} else {	// drop this read as its length is not enough
				++ kstat->dropped[tn];
				continue;
			}
		} else {	// seed not found, now check the tail 2 or 1, if perfect match, drop these 2
			i = wkr->size - 2;
			const char *p = wkr->seq1;
			const char *q = wkr->seq2;
			if( p[i]==kp.adapter_r1[0] && p[i+1]==kp.adapter_r1[1] &&
					q[i]==kp.adapter_r2[0] && q[i+1]==kp.adapter_r2[1] ) {	// a possible hit
				// if it is a real adapter, then Read1 and Read2 should be complimentary
				// in real data, the heading 5 bp are of poor quality, therefoe we test the 6th, 7th
				if( is_revcomp(p[5], q[i-6]) && is_revcomp(q[5], p[i-6]) ) {
					++ kstat->tail_adapter[tn];
					if( i < kp.min_length ) {
						++ kstat->dropped[tn];
						continue;
					}
					CPEREAD_resize( wkr, i );
				}
			} else {	// tail 2 is not good, check tail 1
				++ i;
				if( p[i]==kp.adapter_r1[0] && q[i]==kp.adapter_r2[0] ) {
					if( is_revcomp(p[5], q[i-6]) && is_revcomp(q[5], p[i-6]) && 
						is_revcomp(p[6], q[i-7]) && is_revcomp(q[6], p[i-7]) ) {
						++ kstat->tail_adapter[tn];
						if( i < kp.min_length ) {
							++ kstat->dropped[tn];
							continue;
						}
						CPEREAD_resize( wkr, i );
					}
				}
			}
		}
		writebuffer->b1stored[tn] += sprintf( writebuffer->buffer1[tn]+writebuffer->b1stored[tn],
											"%s%s\n+\n%s\n", wkr->id1, wkr->seq1, wkr->qual1 );
		writebuffer->b2stored[tn] += sprintf( writebuffer->buffer2[tn]+writebuffer->b2stored[tn],
											"%s%s\n+\n%s\n", wkr->id2, wkr->seq2, wkr->qual2 );
	}

	delete [] seed;
}


// producer
int producer_pe_task(const ktrim_param &kp, rabbit::fq::FastqDataPool * fastqPool, FqDataPairChunkQueue& dq){
	vector<string> R1s, R2s;
	extractFileNames( kp.FASTQ1, R1s );
	extractFileNames( kp.FASTQ2, R2s );
	if( R1s.size() != R2s.size() ) {
		fprintf( stderr, "\033[1;31mError: Read1 and Read2 do not contain equal sized files!\033[0m\n" );
		exit(0);
	}
	unsigned int totalFileCount = R1s.size();
	int n_chunks = 0;
	for(unsigned int fileCnt = 0; fileCnt < totalFileCount; fileCnt++){
		rabbit::fq::FastqFileReader * fqFileReader;
		// 判断是否是gz文件 假设双端数据的两个输入文件的类型一致
		bool file_is_gz = false;
		register unsigned int i = R1s[fileCnt].size() - 3;
		register const char * p = R1s[fileCnt].c_str();
		if( p[i]=='.' && p[i+1]=='g' && p[i+2]=='z' ) {
			file_is_gz = true;
			fqFileReader = new rabbit::fq::FastqFileReader(R1s[fileCnt],*fastqPool,R2s[fileCnt],true);
		} else {
			fqFileReader = new rabbit::fq::FastqFileReader(R1s[fileCnt],*fastqPool,R2s[fileCnt],false);
		}
		while (true){
			rabbit::fq::FastqPairChunk *fqPairChunk = new rabbit::fq::FastqPairChunk;
			fqPairChunk->chunk = fqFileReader->readNextPairChunk();
			if(fqPairChunk->chunk == NULL) break;
			dq.Push(n_chunks,fqPairChunk->chunk);
			n_chunks++;
		}
		delete fqFileReader;
	}
	dq.SetCompleted();
	// std::cout<<"Input files have "<<n_chunks<<" chunks "<<std::endl;
	return 0;

}

// comusmer task
void consumer_pe_task(ktrim_param &kp, rabbit::fq::FastqDataPool *fastqPool, FqDataPairChunkQueue &dq, int thread_num, ktrim_stat *kstat,
						unsigned int* buffer_size, writeBufferTotal * buffer_head){
	rabbit::int64 read_count  = 0;
	rabbit::int64 chunk_id;
	rabbit::fq::FastqPairChunk* fqPairChunk = new rabbit::fq::FastqPairChunk;
	CPEREAD * read = new CPEREAD [ READS_PER_BATCH ];
	char * read_data = new char [ MEM_PE_READSET ];
	for(int i = 0, j= 0; i < READS_PER_BATCH ; i++){
		read[i].id1   = read_data + j;
		j += MAX_READ_ID;
		read[i].seq1  = read_data + j;
		j += MAX_READ_CYCLE;
		read[i].qual1 = read_data + j;
		j += MAX_READ_CYCLE;
		
		read[i].id2   = read_data + j;
		j += MAX_READ_ID;
		read[i].seq2  = read_data + j;
		j += MAX_READ_CYCLE;
		read[i].qual2 = read_data + j;
		j += MAX_READ_CYCLE;
	}
	while(dq.Pop(chunk_id,fqPairChunk->chunk)){
		chunkFormat_PE(fqPairChunk->chunk,read,true);
		//找到合适的buffer的位置
		lock(buffer_size_mtx);
		workingThread_PE_C( thread_num, 0, end, read, kstat, &writebuffer, kp );
	}

}

int process_PE_C( const ktrim_param &kp ) {
	rabbit::fq::FastqDataPool * fastqPool = new rabbit::fq::FastqDataPool(256,MEM_PER_CHUNK);
	FqDataPairChunkQueue queue1(256,1);

	unsigned int consumer_num = kp.thread;

	ktrim_stat kstat;
	kstat.dropped	   = new unsigned int [ consumer_num ];
	kstat.real_adapter = new unsigned int [ consumer_num ];
	kstat.tail_adapter = new unsigned int [ consumer_num ];

	// 初始化输出
	unsigned int buffer_size = CHUNK_NUM;
	writeBufferTotal * buffer_head = new writeBufferTotal(CHUNK_NUM,MEM_PER_CHUNK,true);

	// producer
	std::thread producer(producer_pe_task,kp,fastqPool,std::ref(queue1));

	// consumer
	std::thread **consumer_threads = new std::thread* [ consumer_num ]; 
	for(int tn = 0; tn < consumer_num; tn++){
		kstat.dropped[tn] = 0;
		kstat.real_adapter[tn] = 0;
		kstat.tail_adapter[tn] = 0;

		consumer_threads[tn] = new std::thread(std::bind(&consumer_pe_task,kp,fastqPool,std::ref(queue1), tn, &kstat,buffer_size,buffer_head));
	}

	


	string fileName = kp.outpre;
	fileName += ".read1.fq";
	FILE *fout1 = fopen( fileName.c_str(), "wt" );
	fileName[ fileName.size()-4 ] = '2';	// read1 -> read2
	FILE *fout2 = fopen( fileName.c_str(), "wt" );
	if( fout1==NULL || fout2==NULL ) {
		fprintf( stderr, "\033[1;31mError: write file failed!\033[0m\n" );
		fclose( fout1 );
		fclose( fout2 );
		return 103;
	}

	register unsigned int line = 0;
	

	// write trim.log
	fileName = kp.outpre;
	fileName += ".trim.log";
	ofstream fout( fileName.c_str() );
	if( fout.fail() ) { 
		fprintf( stderr, "\033[1;34mError: cannot write log file!\033[0m\n" );
		return 105;
	}
	int dropped_all=0, real_all=0, tail_all=0;
	for( unsigned int i=0; i!=kp.thread; ++i ) {
		dropped_all += kstat.dropped[i];
		real_all += kstat.real_adapter[i];
		tail_all += kstat.tail_adapter[i];
	}
	fout << "Total\t"    << line		<< '\n'
		 << "Dropped\t"  << dropped_all << '\n'
		 << "Aadaptor\t" << real_all	<< '\n'
		 << "TailHit\t"  << tail_all	<< '\n';
	fout.close();

	//free memory
	delete [] kstat.dropped;
	delete [] kstat.real_adapter;
	delete [] kstat.tail_adapter;


	return 0;
}

int process_single_thread_PE_C( const ktrim_param &kp ) {
//	fprintf( stderr, "process_single_thread_PE_C\n" );
	// IO speed-up
	ios::sync_with_stdio( false ); // 关闭cin/cout 与 stdio的同步 提升cin、cout的效率
//	cin.tie( NULL );

	CPEREAD *read = new CPEREAD[ READS_PER_BATCH ];
	register char *read_data = new char[ MEM_PE_READSET ];
	
	for( register int i=0, j=0; i!=READS_PER_BATCH; ++i ) {
		read[i].id1   = read_data + j;
		j += MAX_READ_ID;
		read[i].seq1  = read_data + j;
		j += MAX_READ_CYCLE;
		read[i].qual1 = read_data + j;
		j += MAX_READ_CYCLE;
		
		read[i].id2   = read_data + j;
		j += MAX_READ_ID;
		read[i].seq2  = read_data + j;
		j += MAX_READ_CYCLE;
		read[i].qual2 = read_data + j;
		j += MAX_READ_CYCLE;
	}

	ktrim_stat kstat;
	kstat.dropped	   = new unsigned int [ 1 ];
	kstat.real_adapter = new unsigned int [ 1 ];
	kstat.tail_adapter = new unsigned int [ 1 ];
	kstat.dropped[0] = 0;
	kstat.real_adapter[0] = 0;
	kstat.tail_adapter[0] = 0;

	// buffer for storing the modified reads per thread
	writeBuffer writebuffer;
	writebuffer.buffer1  = new char * [ 1 ];
	writebuffer.buffer2  = new char * [ 1 ];
	writebuffer.b1stored = new unsigned int	[ 1 ];
	writebuffer.b2stored = new unsigned int	[ 1 ];
	writebuffer.buffer1[0] = new char[ BUFFER_SIZE_PER_BATCH_READ_ST ];
	writebuffer.buffer2[0] = new char[ BUFFER_SIZE_PER_BATCH_READ_ST ];

// deal with multiple input files
	vector<string> R1s, R2s;
	extractFileNames( kp.FASTQ1, R1s );
	extractFileNames( kp.FASTQ2, R2s );

	if( R1s.size() != R2s.size() ) {
		fprintf( stderr, "\033[1;31mError: Read1 and Read2 do not contain equal sized files!\033[0m\n" );
		return 110;
	}
	unsigned int totalFiles = R1s.size();
	//cout << "\033[1;34mINFO: " << totalFiles << " paired fastq files will be loaded.\033[0m\n";

	FILE *fout1, *fout2;
	string fileName = kp.outpre;
	fileName += ".read1.fq";
	fout1 = fopen( fileName.c_str(), "wt" );
	fileName[ fileName.size()-4 ] = '2';	// read1 -> read2
	fout2 = fopen( fileName.c_str(), "wt" );
	if( fout1==NULL || fout2==NULL ) {
		fprintf( stderr, "\033[1;31mError: write file failed!\033[0m\n" );
		fclose( fout1 );
		fclose( fout2 );
		return 103;
	}

	register unsigned int line = 0;
	for( unsigned int fileCnt=0; fileCnt!=totalFiles; ++ fileCnt ) {
		bool file_is_gz = false;
		FILE *fq1, *fq2;
		gzFile gfp1, gfp2;
		register unsigned int i = R1s[fileCnt].size() - 3;
		register const char * p = R1s[fileCnt].c_str();
		register const char * q = R2s[fileCnt].c_str();
		if( p[i]=='.' && p[i+1]=='g' && p[i+2]=='z' ) {
			file_is_gz = true;
			gfp1 = gzopen( p, "r" );
			gfp2 = gzopen( q, "r" );
			if( gfp1==NULL || gfp2==NULL ) {
				fprintf( stderr, "\033[1;31mError: open fastq file failed!\033[0m\n" );
				fclose( fout1 );
				fclose( fout2 );
				return 104;
			}
		} else {
			fq1 = fopen( p, "rt" );
			fq2 = fopen( q, "rt" );
			if( fq1==NULL || fq2==NULL ) {
				fprintf( stderr, "\033[1;31mError: open fastq file failed!\033[0m\n" );
				fclose( fout1 );
				fclose( fout2 );
				return 104;
			}
		}

		register unsigned int last_seed;
		vector<unsigned int> seed;
		vector<unsigned int> :: iterator it;

		while( true ) {
			// get fastq reads
			unsigned int loaded;

			if( file_is_gz ) {
				loaded = load_batch_data_PE_GZ( gfp1, gfp2, read, READS_PER_BATCH_ST );
			} else {
				loaded = load_batch_data_PE_C( fq1, fq2, read, READS_PER_BATCH_ST );
			}

			if( loaded == 0 )
				break;
			
			workingThread_PE_C( 0, 0, loaded, read, &kstat, &writebuffer, kp );

			// write output and update fastq statistics
			fwrite( writebuffer.buffer1[0], sizeof(char), writebuffer.b1stored[0], fout1 );
			fwrite( writebuffer.buffer2[0], sizeof(char), writebuffer.b2stored[0], fout2 );

			line += loaded;
			//cerr << '\r' << line << " reads loaded";

			if( file_is_gz ) {
				if( gzeof( gfp1 ) ) break;
			} else {
				if( feof( fq2 ) ) break;
			}
		}

		if( file_is_gz ) {
			gzclose( gfp1 );
			gzclose( gfp2 );
		} else {
			fclose( fq1 );
			fclose( fq2 );
		}
	}
	fclose( fout1 );
	fclose( fout2 );
	//cerr << "\rDone: " << line << " lines processed.\n";

	// write trim.log
	fileName = kp.outpre;
	fileName += ".trim.log";
	ofstream fout( fileName.c_str() );
	if( fout.fail() ) {
		fprintf( stderr, "\033[1;34mError: cannot write log file!\033[0m\n" );
		return 105;
	}

	fout << "Total\t"	 << line					<< '\n'
		 << "Dropped\t"  << kstat.dropped[0]		<< '\n'
		 << "Aadaptor\t" << kstat.real_adapter[0]	<< '\n'
		 << "TailHit\t"  << kstat.tail_adapter[0]	<< '\n';
	fout.close();

	delete [] read;
	delete [] read_data;

	return 0;
}

