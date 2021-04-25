/**
 * Author: Kun Sun (sunkun@szbl.ac.cn)
 * Date:   Feb, 2020
 * This program is part of the Ktrim package
**/

#include <fstream>
#include <sstream>
#include <algorithm>
#include <thread>
#include <stdlib.h>
#include <memory.h>
#include <omp.h>
#include <zlib.h>
#include "common.h"
#include "util.h"
#include <string>
#include <iostream>
#include <thread>
#include "io/FastxChunk.h"
#include "io/FastxStream.h"
#include "io/DataQueue.h"
#include "io/Formater.h"
#include <functional>
#include <atomic>
#include <mutex>
#include "io/Globals.h"
#include "find_seed.h"
#include <cstdint>
// #include <boost/bind.hpp>

using namespace std;
// 同步输出缓存链表的头节点指针
std::mutex mtx_buffer;
std::mutex  mtx_queue;
int buffer_head_pos = 0;

inline int get_num(writeBufferTotal* cb){
	int count = 0;
	writeBufferTotal* tmp = cb;
	while(tmp){
		cout << "in while" << endl;
		count ++; 
		tmp = tmp->next1;
	}
	return count;
}


void inline CSEREAD_resize( CSEREAD * cr, int n ) {
	cr->seq[n] = 0;
	cr->qual[n] = 0;
	cr->size = n;
}

void find_seed( vector<unsigned int> &seed, CSEREAD *read, const ktrim_param & kp ) {
	seed.clear();
	// SE 寻找满足index1的seeds
	register char *poffset  = read->seq;
	register char *indexloc = poffset;
	while( true ) {
		indexloc = strstr( indexloc, kp.adapter_index1 ); // strstr(char*a, char*b)是在a中找到第一次出现b的位置 返回值是char* 找不到时返回NULL
		if( indexloc == NULL )
			break;
		seed.push_back( indexloc - poffset ); // 
		indexloc ++;
	}
	// SE 寻找满足index3的seeds
	poffset  = read->seq + OFFSET_INDEX3;
	indexloc = poffset;
	while( true ) {
		indexloc = strstr( indexloc, kp.adapter_index3 );
		if( indexloc == NULL )
			break;
		seed.push_back( indexloc - poffset ); // 注意：因为这是SE的adapter的第二个index（4-6bp）seed中存储的应该是adapter的第一个bp的位置
		indexloc ++;
	}
	sort( seed.begin(), seed.end());
}

void workingThread_SE_C( unsigned int tn, unsigned int start, unsigned int end, CSEREAD *workingReads,
							ktrim_stat * kstat, writeBuffer * writebuffer, const ktrim_param & kp ) {

//	fprintf( stderr, "=== working thread %d: %d - %d\n", tn, start, end ), "\n";
	// 每次输出缓存区域都将被覆盖
	writebuffer->b1stored[tn] = 0;

	register int i, j;
	register unsigned int last_seed;
	vector<unsigned int> seed;
	vector<unsigned int> :: iterator it;
	const char *p, *q;

	register CSEREAD * wkr = workingReads + start; // wkr 开始指向当前线程第一个read的首地址，表示正在处理中的read
	// 处理当前线程的每一个read
	for( unsigned int ii=start; ii!=end; ++ii, ++wkr ) {
//		fprintf( stderr, "working: %d, %s\n", ii, wkr->id );
		// quality control
		q = wkr->qual;
		j = wkr->size;

		// update in v1.2: support window check
		i = get_quality_trim_cycle_se( q, j, kp );

		if( i == 0 ) { // not long enough
			++ kstat->dropped[ tn ];
			continue;
		}
		if( i != j ) {  //如果进行了quality-trim需要修改CSEREAD的size
			CSEREAD_resize( wkr, i);
		}

		// looking for seed target, 1 mismatch is allowed for these 2 seeds
		// which means seq1 and seq2 at least should take 1 perfect seed match
		find_seed( seed, wkr, kp );
		
		// 处理每个read的所有seed
		last_seed = impossible_seed;	// a position which cannot be in seed 
		for( it=seed.begin(); it!=seed.end(); ++it ) {
			if( *it != last_seed ) {
//				fprintf( stderr, " check seed: %d\n", *it );
			// as there maybe the same value in seq1_seed and seq2_seed,
			// use this to avoid re-calculate that pos
			// 因为SE使用了双索引 并且seed都表示的都是adapter的第一个碱基的位置 所以可能会有重复的seed position 需要跳过
				if( check_mismatch_dynamic_SE_C( wkr, *it, kp ) )
					break;

				last_seed = *it;
			}
		}
		// 退出上面for循环 it!=end时 it表示adapter的位置
		if( it != seed.end() ) {	// adapter found
			++ kstat->real_adapter[tn];
			if( *it >= kp.min_length )	{
				CSEREAD_resize( wkr, *it );
			} else {	// drop this read as its length is not enough
				++ kstat->dropped[tn];
				continue;
			}
		} else {
				// seed not found, now check the tail 2, if perfect match, drop these 2;
				// Single-end reads do not check tail 1
				// 如果所有的seed都不是adapter 则检查最后两个碱基 （因为双索引的长度都是3）判断是否匹配adapter
			i = wkr->size - 2;
			p = wkr->seq;
			if( p[i]==kp.adapter_r1[0] && p[i+1]==kp.adapter_r1[1] ) {
				++ kstat->tail_adapter[tn];
				if( i < kp.min_length ) {
					++ kstat->dropped[tn];
					continue;
				}
				CSEREAD_resize( wkr, i );
			}
		}
		writebuffer->b1stored[tn] += sprintf( writebuffer->buffer1[tn]+writebuffer->b1stored[tn],
												"%s\n%s\n+\n%s\n", wkr->id, wkr->seq, wkr->qual);

	}
}

void workingThread_SE_C_Multi( unsigned int tn, int chunk_id, unsigned int start, unsigned int end, CSEREAD *workingReads,
							ktrim_stat * kstat, writeBufferTotal * writebuffer_total, const ktrim_param & kp ) {
	
	char index_10 = kp.adapter_index1[0];
	char index_11 = kp.adapter_index1[1];
	char index_12 = kp.adapter_index1[2];
	char index_20 = kp.adapter_index1[0];
	char index_21 = kp.adapter_index1[1];
	char index_22 = kp.adapter_index1[2];

	const int16_t ** seed_pos_1 = new const int16_t*[72]; // 9*8
	const int16_t ** seed_pos_2 = new const int16_t*[72]; // 9*8
	for(int i = 0; i < 72;i++){
		seed_pos_1[i] = NULL;
		seed_pos_2[i] = NULL;
	}

//	fprintf( stderr, "=== working thread %d: %d - %d\n", tn, start, end ), "\n";

	// 记录writebuffer_total[chunk_id]写入的位置
	int b1stored = 0;

	register int i, j;
	register unsigned int last_seed;
	// vector<unsigned int> seed;
	// seed.reserve(30);
	// vector<unsigned int> :: iterator it;
	const char *p, *q;

	register CSEREAD * wkr = workingReads + start; // wkr 开始指向当前线程第一个read的首地址，表示正在处理中的read
	// 处理当前线程的每一个read
	for( unsigned int ii=start; ii!=end; ++ii, ++wkr ) {
//		fprintf( stderr, "working: %d, %s\n", ii, wkr->id );
		// quality control
		q = wkr->qual;
		j = wkr->size;

		// update in v1.2: support window check
		i = get_quality_trim_cycle_se( q, j, kp );

		if( i == 0 ) { // not long enough
			++ kstat->dropped[ tn ];
			continue;
		}
		if( i != j ) {  //如果进行了quality-trim需要修改CSEREAD的size
			CSEREAD_resize( wkr, i);
		}

		// looking for seed target, 1 mismatch is allowed for these 2 seeds
		// which means seq1 and seq2 at least should take 1 perfect seed match
		// find_seed( seed, wkr, kp );
		int num_ = wkr->size / 62;
		if(wkr->size % 62 > 2) num_++; 
		find_seed(wkr->seq,num_,seed_pos_1, seed_pos_2, 0, 0, 
					index_10, index_11, index_12, index_20, index_21, index_22);
		
		last_seed = impossible_seed;	// a position which cannot be in seed
		bool isFind = false;
		int fs = 0;
		int16_t tmp;
		for(fs = 0; fs < (num_ << 3); fs++){
			int num_1 = seed_pos_1[fs][0] + 1;
			for(int ss = 1; ss < num_1; ss++){
				tmp = seed_pos_1[fs][ss];
				if(tmp == last_seed || tmp < 0 || tmp >= wkr->size) continue;
				last_seed = tmp;
				if( check_mismatch_dynamic_SE_C( wkr, tmp, kp ) ){
					isFind = true;
					break;
				}
			}
			if(isFind) break;
			int num_2 = seed_pos_2[fs][0] + 1;
			for(int ss = 1; ss < num_2; ss++){
				tmp = seed_pos_2[fs][ss];
				if(tmp == last_seed || tmp < 0 || tmp >= wkr->size) continue;
				last_seed = tmp;
				if( check_mismatch_dynamic_SE_C( wkr, tmp, kp ) ){
					isFind = true;
					break;
				}
			}
			if(isFind) break;
		}

		if(fs < (num_ << 3)){
			++ kstat->real_adapter[tn];
			if( tmp >= kp.min_length )	{
				CSEREAD_resize( wkr, tmp );
			} else {	// drop this read as its length is not enough
				++ kstat->dropped[tn];
				continue;
			}
		}else{
			// 如果所有的seed都不是adapter 则检查最后两个碱基 （因为双索引的长度都是3）判断是否匹配adapter
			i = wkr->size - 2;
			p = wkr->seq;
			if( p[i]==index_10 && p[i+1]==index_11 ) {
				++ kstat->tail_adapter[tn];
				if( i < kp.min_length ) {
					++ kstat->dropped[tn];
					continue;
				}
				CSEREAD_resize( wkr, i );
			}
		}
		
		// delete []seed_pos_1;
		// delete []seed_pos_2;
		int lname = strlen(wkr->id);
		strcpy(writebuffer_total->buffer1[chunk_id]+b1stored,wkr->id);
		b1stored+= lname;
		writebuffer_total->buffer1[chunk_id][b1stored] = '\n';
		b1stored++;
		strcpy(writebuffer_total->buffer1[chunk_id]+b1stored,wkr->seq);
		b1stored+= wkr->size;
		writebuffer_total->buffer1[chunk_id][b1stored] = '\n';
		b1stored++;
		strcpy(writebuffer_total->buffer1[chunk_id]+b1stored,"+\n");
		b1stored += 2;
		strcpy(writebuffer_total->buffer1[chunk_id]+b1stored,wkr->qual);
		b1stored+= wkr->size;
		writebuffer_total->buffer1[chunk_id][b1stored] = '\n';
		b1stored++;
		// b1stored += sprintf( writebuffer->buffer1[tn]+writebuffer->b1stored[tn],
		// 										"%s\n%s\n+\n%s\n", wkr->id, wkr->seq, wkr->qual);
	}
	writebuffer_total->buffer1[chunk_id][b1stored] = 0;
	writebuffer_total->buffer1[chunk_id][MEM_PER_CHUNK-1] = 'd'; // 最后一位设置为标志位 表示当前chunk处理完成
	delete []seed_pos_1;
	delete []seed_pos_2;
}

void write_buffer_producer_se_task(WriteBufferDataPool *bufferPool,WriteBufferQueue &bufferQueue,std::atomic_bool *consumer_task_finished){
	rabbit::int64 id = 0;
	while (!(*consumer_task_finished)){
		writeBufferTotal* buffer = (writeBufferTotal*)malloc(sizeof(writeBufferTotal*));
		bufferPool->Acquire(buffer);
		bufferQueue.Push(id,buffer);
		id++;
	}
	bufferQueue.SetCompleted();
	// 释放buffer_queue没有用到的 buffer
	rabbit::int64 i;
	writeBufferTotal* buffer_ = new writeBufferTotal;
	lock_guard<std::mutex> lg(mtx_queue);
	while (!bufferQueue.IsEmpty()){
		bufferQueue.Pop(i,buffer_);
		bufferPool->Release(buffer_);
	}
}

int producer_se_task(const ktrim_param &kp,rabbit::fq::FastqDataPool* fastqPool,FqDataChunkQueue& dq){
	vector<string> R1s;
	extractFileNames(kp.FASTQU,R1s);
	unsigned int totalFileCount = R1s.size();
	int n_chunks = 0;
	for(unsigned int fileCnt = 0; fileCnt < totalFileCount; fileCnt++){
		rabbit::fq::FastqFileReader * fqFileReader;
		// 判断是否是gz文件
		bool file_is_gz = false;
		register unsigned int i = R1s[fileCnt].size() - 3;
		register const char * p = R1s[fileCnt].c_str();
		if( p[i]=='.' && p[i+1]=='g' && p[i+2]=='z' ) {
			file_is_gz = true;
			fqFileReader = new rabbit::fq::FastqFileReader(R1s[fileCnt],*fastqPool,"",true);
		} else {
			fqFileReader = new rabbit::fq::FastqFileReader(R1s[fileCnt],*fastqPool,"",false);
		}
		while (true){
			rabbit::fq::FastqChunk *fqChunk = new rabbit::fq::FastqChunk;
			fqChunk->chunk = fqFileReader->readNextChunk();
			if(fqChunk->chunk == NULL) break;
			dq.Push(n_chunks,fqChunk->chunk);
			n_chunks++;
		}
		delete fqFileReader;
	}
	dq.SetCompleted();
	// std::cout<<"Input files have "<<n_chunks<<" chunks "<<std::endl;
	return 0;
}



void consumer_se_task(const ktrim_param &kp,rabbit::fq::FastqDataPool * fastqPool, FqDataChunkQueue& dq,
						int th_n, ktrim_stat* kstat,
						std::atomic_int *write_size, std::mutex *mtx,
						std::atomic_int *consumer_finished_num,int consumer_num, std::atomic_bool *consumer_task_finished,
						rabbit::int64 * read_count_total, WriteBufferQueue & bufferQueue){
	
	CSEREAD *readA = new CSEREAD[ READS_PER_BATCH ];
	register char *readA_data = new char[ MEM_SE_READSET ];
	for( register int i=0, j=0; i!=READS_PER_BATCH; ++i ) {
		readA[i].id   = readA_data + j;
		j += MAX_READ_ID;
		readA[i].seq  = readA_data + j;
		j += MAX_READ_CYCLE;
		readA[i].qual = readA_data + j;
		j += MAX_READ_CYCLE;
	}
	

	rabbit::int64 read_count_th = 0; //线程th读取的read总数
	rabbit::int64 id = 0; // chunk_id
	rabbit::int64 read_count = 0; // 当前chunk中的read的数量
	rabbit::fq::FastqChunk* fqChunk  = new rabbit::fq::FastqChunk;
	writeBufferTotal* cur_buffer = buffer_head;
	writeBufferTotal *writebuffer_total_new = new writeBufferTotal;
	while (dq.Pop(id,fqChunk->chunk)){
		read_count = chunkFormat(fqChunk->chunk,readA,true);
		read_count_th += read_count;
		// std::cout<<"the chunk has "<< read_count << " read"<<std::endl;
		fastqPool->Release(fqChunk->chunk);	
		// 判断chunk_id 与 write_size 关系
		(*mtx).lock();
		if(id >= *write_size){
			// 扩容
			// writeBufferTotal *writebuffer_total_new  = new writeBufferTotal(CHUNK_NUM,MEM_PER_CHUNK);
			rabbit::int64 buffer_id;
			bufferQueue.Pop(buffer_id,writebuffer_total_new);

			lock_guard<std::mutex> lg(mtx_buffer);
			cur_buffer = buffer_head;
			while (cur_buffer->next1){
				cur_buffer = cur_buffer ->next1;
			}
			cur_buffer->next1 = writebuffer_total_new;
			cur_buffer = cur_buffer->next1;
			*write_size += writebuffer_total_new->chunk_num_;
		}else{
			// 确定向哪个writeBufferTotal中写入
			int buffer_pos = id / CHUNK_NUM;
			lock_guard<std::mutex> lg(mtx_buffer);
			cur_buffer = buffer_head;
			for (int i = buffer_head_pos; i < buffer_pos; i++){
				cur_buffer = cur_buffer->next1;
			}
		}
		(*mtx).unlock();
		workingThread_SE_C_Multi(th_n, id % cur_buffer->chunk_num_ , 0, read_count, readA, kstat,cur_buffer, kp );
	}
	read_count_total[th_n] = read_count_th;
	// free memory
	delete [] readA;
	delete [] readA_data;
	// cout<<"thread "<<th_n<<" of consumer is done"<<endl; 
	(*consumer_finished_num)++;
	if(*consumer_finished_num == consumer_num){
		*consumer_task_finished = true;
	}
		
}


void writer_se_task(std::atomic_bool *consumer_task_finished, const ktrim_param& kp,
					WriteBufferDataPool* bufferPool, WriteBufferQueue& bufferQueue){
	
	FILE *fout1;
	string fileName = kp.outpre;
	fileName += ".read1.fq";
	fout1 = fopen( fileName.c_str(), "wt" );
	if( fout1==NULL ) {
		fprintf( stderr, "\033[1;31mError: write file failed!\033[0m\n" );
		fclose( fout1 );
		return;
	}
	
	writeBufferTotal * cur_buffer = buffer_head;
	int cur_pos = 0;
	int write_pos = 0;  // chunk_id % cur_buffer->chunk_num_
	int chunk_num = cur_buffer->chunk_num_;
	
	while(!(*consumer_task_finished)){
		while(cur_buffer->buffer1[write_pos][MEM_PER_CHUNK-1] != 'd'){
			if(cur_buffer->buffer1[write_pos][MEM_PER_CHUNK-1] == 'd' || (*consumer_task_finished))
				break;
		}
		if(!(*consumer_task_finished)){
			fprintf(fout1,cur_buffer->buffer1[write_pos]);
			write_pos++;

			if(write_pos >= chunk_num){
				write_pos = write_pos % chunk_num; // 0 
				while(cur_buffer->next1 == NULL){
					if(cur_buffer->next1 != NULL || (*consumer_task_finished)){
						break;
					}
				}
				if (*consumer_task_finished){
					cur_buffer = cur_buffer->next1;
					break;
				}
				lock_guard<std::mutex> lg(mtx_buffer);
				writeBufferTotal * tmp = cur_buffer;
				cur_buffer = cur_buffer->next1;
				chunk_num  = cur_buffer->chunk_num_;
				buffer_head  = cur_buffer;
				buffer_head_pos++;
				bufferPool->Release(tmp);
			} // change next buffer over
		}
	}

	// 释放bufferQueue中用不到的buffer
	lock_guard<std::mutex> lg(mtx_queue);
	rabbit::int64 int_a;
	writeBufferTotal * buffer_a = new writeBufferTotal;
	while (!bufferQueue.IsEmpty()){
		bufferQueue.Pop(int_a,buffer_a);
		bufferPool->Release(buffer_a);
	}
	

	// if(read_pos >= CHUNK_NUM){
	// 	read_pos = read_pos % CHUNK_NUM;
	// 	if (cur_buffer->next1 == NULL) {
	// 		fclose(fout1);
	// 		return;
	// 	}	
	// 	cur_buffer = cur_buffer->next1;
	// }
	if (cur_buffer == NULL){
		fclose(fout1);
		return;
	}

	while (cur_buffer->buffer1[write_pos][MEM_PER_CHUNK-1] == 'd'){
		fprintf(fout1,cur_buffer->buffer1[write_pos]);
		write_pos++;
		if(write_pos >= chunk_num){
			write_pos = write_pos % chunk_num; // = 0 
			if (cur_buffer->next1 == NULL) return;
			writeBufferTotal * tmp = cur_buffer;
			cur_buffer = cur_buffer->next1;
			chunk_num = cur_buffer->chunk_num_;
			lock_guard<std::mutex> lg(mtx_buffer);
			buffer_head_pos++;
			buffer_head = cur_buffer;
			bufferPool->Release(tmp);
		}
	}
	fclose(fout1);
}

int process_SE_C( const ktrim_param &kp ) {

	rabbit::fq::FastqDataPool * fastqPool = new rabbit::fq::FastqDataPool(256,MEM_PER_CHUNK);
	FqDataChunkQueue queue1(256,1);
	//std::atomic_int last_chunk_id(-1); // 保证输入和输出的有序性 k号线程只能处理last_chunk_id % consumer_num = k 的chunk
	unsigned int consumer_num = kp.thread;
	std::atomic_int consumer_finished_num;
	std::atomic_bool consumer_task_finished;
	std::atomic_init(&consumer_finished_num,0);
	std::atomic_init(&consumer_task_finished,false);

	ktrim_stat kstat;
	kstat.dropped	   = new unsigned int [ consumer_num ];
	kstat.real_adapter = new unsigned int [ consumer_num ];
	kstat.tail_adapter = new unsigned int [ consumer_num ];
	rabbit::int64 * read_count_total = new rabbit::int64[consumer_num];

	for(unsigned int i=0; i!=consumer_num; ++i) {
		kstat.dropped[i] = 0;
		kstat.real_adapter[i] = 0;
		kstat.tail_adapter[i] = 0;
	}

	// write buffer producer
	WriteBufferDataPool* bufferPool = new WriteBufferDataPool(2,MEM_PER_CHUNK);
	WriteBufferQueue bufferQueue(2,1);
	std::thread write_buffer_producer(std::bind(&write_buffer_producer_se_task,bufferPool,std::ref(bufferQueue),&consumer_task_finished));
	rabbit::int64 buffer_id;
	bufferQueue.Pop(buffer_id,buffer_head);
	
	std::atomic_int write_size; //记录当前输出buffer的大小
	std::atomic_init(&write_size,buffer_head->chunk_num_);
	std::mutex mtx;

	// producer
	std::thread producer(producer_se_task,kp,fastqPool,std::ref(queue1));
	
	// consumer
	std::thread ** threads = new std::thread* [ consumer_num ];
	for(int t = 0;t < consumer_num; t++){
		threads[t] = new std::thread(std::bind(&consumer_se_task, kp, fastqPool, std::ref(queue1), 
										t, &kstat, &write_size, &mtx, &consumer_finished_num, 
										consumer_num, &consumer_task_finished, read_count_total,
										std::ref(bufferQueue)));
	}

	// writer
	std::thread* writer = NULL;
	writer = new std::thread(std::bind(&writer_se_task,&consumer_task_finished,kp, bufferPool, std::ref(bufferQueue)));
	
	// thread join
	write_buffer_producer.join();
	producer.join();
	for (int t = 0; t < consumer_num; t++){
		threads[t]->join();
	}
	writer->join();

	
	// write trim.log
	string fileName = kp.outpre;
	fileName += ".trim.log";
	ofstream fout( fileName.c_str() );
	if( fout.fail() ) { 
		fprintf( stderr, "\033[1;34mError: cannot write log file!\033[0m\n" );
		return 105;
	}
	int dropped_all=0, real_all=0, tail_all=0;
	rabbit::int64 line = 0;
	for( unsigned int i=0; i!=consumer_num; ++i ) {
		line += read_count_total[i];
		dropped_all += kstat.dropped[i];
		real_all += kstat.real_adapter[i];
		tail_all += kstat.tail_adapter[i];
	}
	fout << "Total\t"	<< line		<< '\n'
		 << "Dropped\t" << dropped_all << '\n'
		 << "Aadaptor\t" << real_all	<< '\n'
		 << "TailHit\t" << tail_all	<< '\n';
	fout.close();

	// free memory
	delete [] kstat.dropped;
	delete [] kstat.real_adapter;
	delete [] kstat.tail_adapter;
	delete [] read_count_total;
	delete bufferPool;

	for(int t = 0; t < consumer_num; t++){
		delete threads[t];
	}
	delete [] threads;
	
	return 0;
}
	
	
	
	
	
	
/* 
int process_multi_thread_SE_C( const ktrim_param &kp ) {
	// in this version, two data containers are used and auto-swapped for working and loading data
	CSEREAD *readA = new CSEREAD[ READS_PER_BATCH ];
	CSEREAD *readB = new CSEREAD[ READS_PER_BATCH ];

	register char *readA_data = new char[ MEM_SE_READSET ];
	register char *readB_data = new char[ MEM_SE_READSET ];

	for( register int i=0, j=0; i!=READS_PER_BATCH; ++i ) {
		readA[i].id   = readA_data + j;
		readB[i].id   = readB_data + j;
		j += MAX_READ_ID;
		readA[i].seq  = readA_data + j;
		readB[i].seq  = readB_data + j;
		j += MAX_READ_CYCLE;
		readA[i].qual = readA_data + j;
		readB[i].qual = readB_data + j;
		j += MAX_READ_CYCLE;
	}

	CSEREAD *workingReads, *loadingReads, *swapReads;

	ktrim_stat kstat;
	kstat.dropped	   = new unsigned int [ kp.thread ];
	kstat.real_adapter = new unsigned int [ kp.thread ];
	kstat.tail_adapter = new unsigned int [ kp.thread ];

	// buffer for storing the modified reads per thread
	writeBuffer writebuffer;
	writebuffer.buffer1  = new char * [ kp.thread ];
	writebuffer.b1stored = new unsigned int	[ kp.thread ];

	for(unsigned int i=0; i!=kp.thread; ++i) {
		writebuffer.buffer1[i] = new char[ BUFFER_SIZE_PER_BATCH_READ ];

		kstat.dropped[i] = 0;
		kstat.real_adapter[i] = 0;
		kstat.tail_adapter[i] = 0;
	}

	// deal with multiple input files
	vector<string> R1s;
	extractFileNames( kp.FASTQU, R1s );
	unsigned int totalFiles = R1s.size();
	//cout << "\033[1;34mINFO: " << totalFiles << " single-end fastq files will be loaded.\033[0m\n";

	FILE *fout1;
	string fileName = kp.outpre;
	fileName += ".read1.fq";
	fout1 = fopen( fileName.c_str(), "wt" );
	if( fout1==NULL ) {
		fprintf( stderr, "\033[1;31mError: write file failed!\033[0m\n" );
		fclose( fout1 );
		return 103;
	}

	register unsigned int line = 0;
	unsigned int threadCNT = kp.thread - 1;
	for( unsigned int fileCnt=0; fileCnt!=totalFiles; ++ fileCnt ) {
		bool file_is_gz = false;
		FILE *fq;
		gzFile gfp;
		register unsigned int i = R1s[fileCnt].size() - 3;
		register const char * p = R1s[fileCnt].c_str();
		if( p[i]=='.' && p[i+1]=='g' && p[i+2]=='z' ) {
			file_is_gz = true;
			gfp = gzopen( p, "r" );
			if( gfp == NULL ) {
				fprintf( stderr, "\033[1;31mError: open fastq file failed!\033[0m\n" );
				fclose( fout1 );
				return 104;
			}
		} else {
			fq = fopen( p, "rt" );
			if( fq == NULL ) {
				fprintf( stderr, "\033[1;31mError: open fastq file failed!\033[0m\n" );
				fclose( fout1 );
				return 104;
			}
		}
		// initialization
		// get first batch of fastq reads
		unsigned int loaded;
		bool metEOF;
		if( file_is_gz ) {
			loaded = load_batch_data_SE_GZ( gfp, readA, READS_PER_BATCH );
			metEOF = gzeof( gfp );
		} else {
			loaded = load_batch_data_SE_C( fq, readA, READS_PER_BATCH );
			metEOF = feof( fq );
		}
		if( loaded == 0 ) break; // TODO 这里是不是应该用continue？
		// if( loaded == 0 ) continue;
//		fprintf( stderr, "Loaded %d, metEOF=%d\n", loaded, metEOF );

		loadingReads = readB;
		workingReads = readA;
		bool nextBatch = true;
		unsigned int threadLoaded;
		unsigned int NumWkThreads; // 数据分析的线程数量
		while( nextBatch ) {
			// start parallalization
			omp_set_num_threads( kp.thread );
			#pragma omp parallel
			{
				unsigned int tn = omp_get_thread_num(); // tn是当前线程的编号
				// if EOF is met, then all threads are used for analysis
				// otherwise 1 thread will do data loading
				if( metEOF ) {
					NumWkThreads = kp.thread;
					unsigned int start = loaded * tn / kp.thread;
					unsigned int end   = loaded * (tn+1) / kp.thread;
					workingThread_SE_C( tn, start, end, workingReads, &kstat, &writebuffer, kp );
					nextBatch = false;
				} else {	// use 1 thread to load file, others for trimming
					NumWkThreads = threadCNT;  //threadCNT = kp.thread -1 ;

					if( tn == threadCNT ) { // 最后一个线程进行data load
						if( file_is_gz ) {
							threadLoaded = load_batch_data_SE_GZ( gfp, loadingReads, READS_PER_BATCH );
							metEOF = gzeof( gfp );
						} else {
							threadLoaded = load_batch_data_SE_C( fq, loadingReads, READS_PER_BATCH );
							metEOF = feof( fq );
						}
						nextBatch = (threadLoaded!=0);
//						fprintf( stderr, "Loaded %d, metEOF=%d\n", threadLoaded, metEOF );
						//cerr << "Loading thread: " << threadLoaded << ", " << metEOF << ", " << nextBatch << '\n';
					} else {
						// 其他线程进行数据处理 
						unsigned int start = loaded * tn / threadCNT;
						unsigned int end   = loaded * (tn+1) / threadCNT;
						workingThread_SE_C( tn, start, end, workingReads, &kstat, &writebuffer, kp );
					}
				}
			} // parallel body
			// swap workingReads and loadingReads for next loop
			swapReads	= loadingReads;
			loadingReads = workingReads;
			workingReads = swapReads;
			// write output and update fastq statistics
			for( unsigned int ii=0; ii!=NumWkThreads; ++ii ) {
				fwrite( writebuffer.buffer1[ii], sizeof(char), writebuffer.b1stored[ii], fout1 );
			}
			line += loaded;
			loaded = threadLoaded;
			//cerr << '\r' << line << " reads loaded";
		}

		if( file_is_gz ) {
			gzclose( gfp );
		} else {
			fclose( fq );
		}
	}

	fclose( fout1 );
	//cerr << "\rDone: " << line << " lines processed.\n";

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
	fout << "Total\t"	<< line		<< '\n'
		 << "Dropped\t" << dropped_all << '\n'
		 << "Aadaptor\t" << real_all	<< '\n'
		 << "TailHit\t" << tail_all	<< '\n';
	fout.close();

	//free memory
	for(unsigned int i=0; i!=kp.thread; ++i) {
		delete writebuffer.buffer1[i];
	}
	delete [] writebuffer.buffer1;
	delete [] kstat.dropped;
	delete [] kstat.real_adapter;
	delete [] kstat.tail_adapter;

	delete [] readA;
	delete [] readB;
	delete [] readA_data;
	delete [] readB_data;

	return 0;
}
*/

int process_single_thread_SE_C( const ktrim_param &kp ) {

	CSEREAD *read = new CSEREAD[ READS_PER_BATCH_ST ];  // 存储每条read各个字段的首地址
	register char *read_data = new char[ MEM_SE_READSET ]; // 存放read数据的内存区域

	// 初始化变量read中的地址
	for( register int i=0, j=0; i!=READS_PER_BATCH; ++i ) {
		read[i].id   = read_data + j;
		j += MAX_READ_ID;
		read[i].seq  = read_data + j;
		j += MAX_READ_CYCLE;
		read[i].qual = read_data + j;
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
	writebuffer.b1stored = new unsigned int	[ 1 ];
	writebuffer.buffer1[0] = new char[ BUFFER_SIZE_PER_BATCH_READ_ST ];

	// deal with multiple input files
	vector<string> R1s; // 所有输入文件的名称
	extractFileNames( kp.FASTQU, R1s );

	unsigned int totalFiles = R1s.size();
	//cout << "\033[1;34mINFO: " << totalFiles << " single-end fastq files will be loaded.\033[0m\n";

	FILE *fout1;
	string fileName = kp.outpre;
	fileName += ".read1.fq";
	fout1 = fopen( fileName.c_str(), "wt" );
	if( fout1==NULL ) {
		fprintf( stderr, "\033[1;31mError: write file failed!\033[0m\n" );
		fclose( fout1 );
		return 103;
	}

	//ifstream fq1; 
	//line 表示已经读取的read的数量
	register unsigned int line = 0;
	for( unsigned int fileCnt=0; fileCnt!=totalFiles; ++ fileCnt ) {
		//fq1.open( R1s[fileCnt].c_str() );
		bool file_is_gz = false;
		FILE *fq;
		gzFile gfp;
		register unsigned int i = R1s[fileCnt].size() - 3;
		register const char * p = R1s[fileCnt].c_str();
		if( p[i]=='.' && p[i+1]=='g' && p[i+2]=='z' ) {
//			fprintf( stderr, "GZ file!\n" );
			file_is_gz = true;
			gfp = gzopen( p, "r" );
			if( gfp == NULL ) {
				fprintf( stderr, "\033[1;31mError: open fastq file failed!\033[0m\n" );
				fclose( fout1 );
				return 104;
			}
		} else {
			fq = fopen( p, "rt" );
			if( fq == NULL ) {
				fprintf( stderr, "\033[1;31mError: open fastq file failed!\033[0m\n" );
				fclose( fout1 );
				return 104;
			}
		}

		// register unsigned int last_seed;
		// vector<unsigned int> seed;
		// vector<unsigned int> :: iterator it;

		// 开始读取每个文件的内容 按照batch读取 batch_zise = READS_PER_BATCH_ST 最多读取batch_zise
		while( true ) {
			// get fastq reads
			//unsigned int loaded = load_batch_data_SE( fq1, read, READS_PER_BATCH_ST );
			unsigned int loaded;
			if( file_is_gz ) {
				// loaded 是 读取的 CSEREAD 的真实数量
				loaded = load_batch_data_SE_GZ( gfp, read, READS_PER_BATCH_ST );
//				fprintf( stderr, "Loaded=%d\n", loaded );
			} else {
				loaded = load_batch_data_SE_C( fq, read, READS_PER_BATCH_ST );
			}

			if( loaded == 0 )
				break;
		
//			fprintf( stderr, "Work\n" );
			workingThread_SE_C( 0, 0, loaded, read, &kstat, &writebuffer, kp );

			// write output and update fastq statistics
//			fprintf( stderr, "Output\n" );
			fwrite( writebuffer.buffer1[0], sizeof(char), writebuffer.b1stored[0], fout1 );

			line += loaded;
			//cerr << '\r' << line << " reads loaded";

			//if( fq1.eof() ) break;
//			fprintf( stderr, "Check\n" );
			if( file_is_gz ) {
				if( gzeof( gfp ) ) break;
			} else {
				if( feof( fq ) ) break;
			}
		}
		//fq1.close();
		if( file_is_gz ) {
			gzclose( gfp );
		} else {
			fclose( fq );
		}
	}
	fclose( fout1 );
	//cerr << "\rDone: " << line << " lines processed.\n";

	// write trim.log
	fileName = kp.outpre;
	fileName += ".trim.log";
	ofstream fout( fileName.c_str() );
	if( fout.fail() ) { 
		fprintf( stderr, "\033[1;34mError: cannot write log file!\033[0m\n" );
		return 105;
	}

	fout << "Total\t"    << line					<< '\n'
		 << "Dropped\t"  << kstat.dropped[0]		<< '\n'
		 << "Aadaptor\t" << kstat.real_adapter[0]	<< '\n'
		 << "TailHit\t"  << kstat.tail_adapter[0]	<< '\n';
	fout.close();

	//free memory
//	delete buffer1;
//	delete buffer2;
//	delete fq1buffer;
	delete [] read;
	delete [] read_data;

	return 0;
}

