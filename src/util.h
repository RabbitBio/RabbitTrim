/**
 * Author: Kun Sun (sunkun@szbl.ac.cn)
 * Date: Feb, 2020
 * This program is part of the Ktrim package
**/

#ifndef _KTRIM_UTIL_
#define _KTRIM_UTIL_

#include <fstream>
#include <sstream>
#include <algorithm>
#include <thread>
#include <ctime>
#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <omp.h>
#include <math.h>
#include <zlib.h>
#include "common.h"
#include<sys/time.h>
#include "io/FastxChunk.h"
#include "io/Reference.h"
#include "io/Formater.h"
#include <stdint.h>
#include <sys/stat.h>
#include <string>
#include <vector>
#include <string.h>
using namespace std;

writeBufferTotal* buffer_head;
writeBufferTotal* buffer_tail;

typedef rabbit::core::TDataQueue<rabbit::fq::FastqDataChunk> FqDataChunkQueue;
typedef rabbit::core::TDataPool<writeBufferTotal> WriteBufferDataPool; //定义write buffer 的数据池
typedef rabbit::core::TDataQueue<writeBufferTotal> WriteBufferQueue;

inline bool startsWith(string const &value, string const &starting) {
  if (starting.size() > value.size()) return false;
  return equal(starting.begin(), starting.end(), value.begin());
}
inline bool endsWith(string const &value, string const &ending) {
  if (ending.size() > value.size()) return false;
  return equal(ending.rbegin(), ending.rend(), value.rbegin());
}

vector<string> split(const string& str, const string& delim) {
	vector<string> res;
	if("" == str) return res;
	//先将要切割的字符串从string类型转换为char*类型
	char * strs = new char[str.length() + 1] ; //不要忘了
	strcpy(strs, str.c_str()); 
 
	char * d = new char[delim.length() + 1];
	strcpy(d, delim.c_str());
 
	char *p = strtok(strs, d);
	while(p) {
		res.emplace_back(p); //存入结果数组
		p = strtok(NULL, d);
	}
 
	return res;
}

string& ClearHeadTailSpace(string &str)   
{  
    if (str.empty())   
    {  
        return str;  
    }  
    str.erase(0,str.find_first_not_of(" "));  
    str.erase(str.find_last_not_of(" ") + 1);  
    return str;  
}  

// Determine if the file exists or not
inline bool isFileExist(const char* filename){
	struct stat buffer;   
  	return (stat(filename, &buffer) == 0); 
}

// format
int neoGetLine(rabbit::fq::FastqDataChunk *&chunk, uint64_t &pos, uint64_t &len) {
  int start_pos = pos;
  const char *data = (char *)chunk->data.Pointer();
  const uint64_t chunk_size = chunk->size + 1;
  /*
  while (pos <= chunk_size && data[pos] != '\n'){
    pos++;
  }
  len = pos - start_pos - 1;
  return len;
  */
  while (pos <= chunk_size) {
    if (data[pos] == '\n') {
      // find a line
      pos++;
      len = pos - start_pos - 1;
      return len;
    } else {
      pos++;
    }
  }
  return 0;
}

int chunkFormat(rabbit::fq::FastqDataChunk *fqDataChunk, CSEREAD *read, bool mHasQuality = true) {
  rabbit::fq::FastqDataChunk *chunk = fqDataChunk;
  uint64_t seq_count = 0;
  uint64_t pos_ = 0;
  neoReference ref;

  while (true) {
    ref.base = chunk->data.Pointer();
    ref.pname = pos_;
    if (neoGetLine(chunk, pos_, ref.lname)) {
      ref.pseq = pos_;
    } else {
      break;
    }
    neoGetLine(chunk, pos_, ref.lseq);
    ref.pstrand = pos_;
    neoGetLine(chunk, pos_, ref.lstrand);
    ref.pqual = pos_;
    neoGetLine(chunk, pos_, ref.lqual);
    
    // print_read(ref);
	// 构造CSEREAD
	memcpy(read[seq_count].id,ref.base+ref.pname,ref.lname);
	memcpy(read[seq_count].seq,ref.base+ref.pseq,ref.lseq);
	memcpy(read[seq_count].qual,ref.base+ref.pqual,ref.lqual);
	(read[seq_count].id)[ref.lname] = 0;
	(read[seq_count].seq)[ref.lseq] = 0;
	(read[seq_count].qual)[ref.lqual] = 0;
	read[seq_count].size = ref.lseq;
	seq_count++;
  }

  return seq_count;
}

unsigned int chunkFormat_PE(rabbit::fq::FastqDataPairChunk* chunk,CPEREAD * read, bool mHasQuality = true){
	rabbit::fq::FastqDataChunk * left_chunk = chunk->left_part;
	rabbit::fq::FastqDataChunk * right_chunk = chunk->right_part;
	unsigned int seq_count = 0;
	uint64_t pos_1 = 0;
	neoReference ref1;
	uint64_t pos_2 = 0;
	neoReference ref2;

	while(true){
		ref1.base = left_chunk->data.Pointer();
		ref2.base = right_chunk->data.Pointer();
		ref1.pname = pos_1;
		ref2.pname = pos_2;
		if (neoGetLine(left_chunk, pos_1, ref1.lname) && neoGetLine(right_chunk, pos_2, ref2.lname)) {
			ref1.pseq = pos_1;
			ref2.pseq = pos_2;
		} else {
			break;
		}
		neoGetLine(left_chunk, pos_1, ref1.lseq);
		neoGetLine(right_chunk,pos_2, ref2.lseq);
		ref1.pstrand = pos_1;
		ref2.pstrand = pos_2;
		neoGetLine(left_chunk, pos_1, ref1.lstrand);
		neoGetLine(right_chunk, pos_2, ref2.lstrand);
		ref1.pqual = pos_1;
		ref2.pqual = pos_2;
		neoGetLine(right_chunk, pos_2, ref2.lqual);
		
		// print_read(ref);
		// 构造CPEREAD
		memcpy(read[seq_count].id1,ref1.base+ref1.pname,ref1.lname);
		memcpy(read[seq_count].seq1,ref1.base+ref1.pseq,ref1.lseq);
		memcpy(read[seq_count].qual1,ref1.base+ref1.pqual,ref1.lqual);
		(read[seq_count].id1)[ref1.lname] = 0;
		(read[seq_count].seq1)[ref1.lseq] = 0;
		(read[seq_count].qual1)[ref1.lqual] = 0;
		read[seq_count].size = ref1.lseq;

		memcpy(read[seq_count].id2,ref2.base+ref2.pname,ref2.lname);
		memcpy(read[seq_count].seq2,ref2.base+ref2.pseq,ref2.lseq);
		memcpy(read[seq_count].qual2,ref2.base+ref2.pqual,ref2.lqual);
		(read[seq_count].id2)[ref2.lname] = 0;
		(read[seq_count].seq2)[ref2.lseq] = 0;
		(read[seq_count].qual2)[ref2.lqual] = 0;
		read[seq_count].size2 = ref2.lseq;
		seq_count++;
	}
	return seq_count;

}


// 多线程时间统计
double getTime(){
	struct timeval tv;
	gettimeofday(&tv,NULL);
	return (double)tv.tv_sec+(double)tv.tv_usec/1000000;
}

// extract file names 
void extractFileNames( const char *str, vector<string> & Rs ) {
	string fileName = "";
	for(unsigned int i=0; str[i]!='\0'; ++i) {
		if( str[i] == FILE_SEPARATOR ) {
			Rs.push_back( fileName );
			fileName.clear();
		} else {
			fileName += str[i];
		}
	}
	if( fileName.size() )   // in case there is a FILE_SEPARATOR at the end
		Rs.push_back( fileName );
}

//load 1 batch of data, using purely C-style
unsigned int load_batch_data_SE_C( FILE * fp, CSEREAD *loadingReads, unsigned int num ) {
	register unsigned int loaded = 0;
//	string unk;
	register CSEREAD *p = loadingReads;
	register CSEREAD *q = p + num;
//	clock_t start = clock();
	while( p != q ) {
		if( fgets( p->id, MAX_READ_ID,  fp ) == NULL ) break; // 当没有下一行时退出
		fgets( p->seq,  MAX_READ_CYCLE, fp );
		fgets( p->qual, MAX_READ_ID,    fp );	// this line is useless
		fgets( p->qual, MAX_READ_CYCLE, fp );

		// remove the tail '\n'
		register unsigned int i = strlen( p->seq ) - 1;
		p->size = i;
		p->seq[i]  = 0;
		p->qual[i] = 0;

		++ p;
	}
//	clock_t end = clock();
//	fprintf( stderr, "%.0f ", (end-start)*1000.0/CLOCKS_PER_SEC );

	return p-loadingReads; // 返回的是数组loadingReads中CSEREAD的数量
}

unsigned int load_batch_data_SE_GZ( gzFile gfp, CSEREAD *loadingReads, unsigned int num ) {
	register unsigned int loaded = 0;
//	string unk;
	register CSEREAD *p = loadingReads;
	register CSEREAD *q = p + num;
	register char c;
	register unsigned int i;

	while( p != q ) {
		if( gzgets( gfp, p->id, MAX_READ_ID ) == NULL ) break; // gzgets函数：从文件中读入MAX_READ_ID个字符或者读入一行 此处是获取某read的ID字段
		gzgets( gfp, p->seq,  MAX_READ_CYCLE );
		gzgets( gfp, p->qual, MAX_READ_ID );	// this line is useless
		gzgets( gfp, p->qual, MAX_READ_CYCLE );
//		fprintf( stderr, "GZload: %s%s%s", p->id, p->seq, p->qual);

		// remove the tail '\n'
		i = strlen( p->seq ) - 1;
		p->size = i;
		p->seq[i]  = 0;
		p->qual[i] = 0;

		++ p;
	}

	return p-loadingReads;
}

unsigned int load_batch_data_PE_C( FILE *fq1, FILE *fq2, CPEREAD *loadingReads, unsigned int num ) {
	//register unsigned int loaded = 0;
	register CPEREAD *p = loadingReads;
	register CPEREAD *q = p + num;
	register CPEREAD *s = p;
	// load read1
//	size_t read_char_num = 0;
//	clock_t start = clock();
	while( p != q ) {
		if( fgets( p->id1, MAX_READ_ID,  fq1 ) == NULL ) break;
		fgets( p->seq1,  MAX_READ_CYCLE, fq1 );
		fgets( p->qual1, MAX_READ_ID,    fq1 );	// this line is useless
		fgets( p->qual1, MAX_READ_CYCLE, fq1 );
		p->size = strlen( p->seq1 );

/*
		register char *pc = p->id1;
		if( getline( &pc, &read_char_num, fq1 ) == -1 ) break;
		pc = p->seq1;
		p->size = getline( &pc, &read_char_num, fq1 );
		pc = p->qual1;
		getline( &pc, &read_char_num, fq1 );
		getline( &pc, &read_char_num, fq1 );
*/

		++ p;
	}
//	clock_t end = clock();
//	fprintf( stderr, "%.1f ", (end-start)*1000.0/CLOCKS_PER_SEC );

	// load read2
//	start = clock();
	while( s != p ) { // 读取和read1中相同数量的数据
		fgets( s->id2,   MAX_READ_ID,    fq2 );
		fgets( s->seq2,  MAX_READ_CYCLE, fq2 );
		fgets( s->qual2, MAX_READ_ID,    fq2 );	// this line is useless
		fgets( s->qual2, MAX_READ_CYCLE, fq2 );
		s->size2 = strlen( s->seq2 );
/*
		register char *pc = s->id2;
		if( getline( &pc, &read_char_num, fq2 ) == -1 ) break;
		pc = s->seq2;
		s->size2 = getline( &pc, &read_char_num, fq2 );
		pc = s->qual2;
		getline( &pc, &read_char_num, fq2 );
		getline( &pc, &read_char_num, fq2 );
*/
		++ s;
	}
//	end = clock();
//	fprintf( stderr, "%.1f ", (end-start)*1000.0/CLOCKS_PER_SEC );

	// update in v1.1.0
	// deal with size matter: if read 1 and read 2 sequences are of different size, keep the shorter one
//	start = clock();
	s = loadingReads;
	while( s != p ) {
		if( s->size > s->size2 )
			s->size = s->size2;

		// remove the tail '\n'
		-- s->size;
		register unsigned int i = s->size;
		s->seq1[  i ] = 0;
		s->qual1[ i ] = 0;
		s->seq2[  i ] = 0;
		s->qual2[ i ] = 0;

		++ s;
	}
//	end = clock();
//	fprintf( stderr, "%.1f\n", (end-start)*1000.0/CLOCKS_PER_SEC );

	return p - loadingReads;
}

unsigned int load_batch_data_PE_GZ( gzFile gfp1, gzFile gfp2, CPEREAD *loadingReads, unsigned int num ) {
	//register unsigned int loaded = 0;
	register CPEREAD *p = loadingReads;
	register CPEREAD *q = p + num;
	register CPEREAD *s = p;
	while( p != q ) {
		if( gzgets( gfp1, p->id1, MAX_READ_ID ) == NULL ) break;
		gzgets( gfp1, p->seq1,  MAX_READ_CYCLE );
		gzgets( gfp1, p->qual1, MAX_READ_CYCLE );	// this line is useless
		gzgets( gfp1, p->qual1, MAX_READ_CYCLE );
		p->size = strlen( p->seq1 );

		++ p;
	}
	while( s != p ) {
		gzgets( gfp2, s->id2,   MAX_READ_ID );
		gzgets( gfp2, s->seq2,  MAX_READ_CYCLE );
		gzgets( gfp2, s->qual2, MAX_READ_CYCLE );	// this line is useless
		gzgets( gfp2, s->qual2, MAX_READ_CYCLE );
		s->size2 = strlen( s->seq2 );

		++ s;
	}

	// update in v1.1.0
	// deal with size matter: if read 1 and read 2 sequences are of different size, keep the shorter one
	s = loadingReads;
	while( s != p ) {
		if( s->size > s->size2 )
			s->size = s->size2;

		// remove the tail '\n'
		-- s->size;
		s->seq1[  s->size] = 0;
		s->qual1[ s->size] = 0;
		s->seq2[  s->size] = 0;
		s->qual2[ s->size] = 0;

		++ s;
	}

	return p - loadingReads;
}

/*
 * use dynamic max_mismatch as the covered size can range from 3 to a large number such as 50,
 * here the maximum mismatch allowed is LEN/8
*/
bool check_mismatch_dynamic_SE_C( CSEREAD *read, unsigned int pos, const ktrim_param & kp ) {
	register unsigned int mis = 0;
	register unsigned int i, len;
	len = read->size - pos;
	if( len > kp.adapter_len )
		len = kp.adapter_len;
	// len =  min(seq_left_len,adapter_len)

	register unsigned int max_mismatch_dynamic; // 允许出现的最大不匹配的碱基数量 
	// update in v1.1.0: allows the users to set the proportion of mismatches
	if( kp.use_default_mismatch ) {
		max_mismatch_dynamic = len >> 3;
		// if( (max_mismatch_dynamic<<3) != len )
		// 	++ max_mismatch_dynamic;
	} else {
		// max_mismatch_dynamic = ceil( len * kp.mismatch_rate );
		max_mismatch_dynamic = floor( len * kp.mismatch_rate );
	}

	register const char *p = read->seq;
	for( i=0; i!=len; ++i ) {
		if( p[pos+i] != kp.adapter_r1[i] ) {
			if( mis == max_mismatch_dynamic ) return false;
			++ mis;
		}
	}

	return true;
}
//bool check_mismatch_dynamic_SE_C( CSEREAD *read, unsigned int pos, const ktrim_param & kp ) {
//	register unsigned int mis = 0;
//	register unsigned int i, len;
//	len = read->size - pos;
//	if( len > kp.adapter_len )
//		len = kp.adapter_len;
//	// len =  min(seq_left_len,adapter_len)
//
//	register unsigned int max_mismatch_dynamic; // 允许出现的最大不匹配的碱基数量 
//	// update in v1.1.0: allows the users to set the proportion of mismatches
//	if( kp.use_default_mismatch ) {
//		max_mismatch_dynamic = len >> 3;
//		if( (max_mismatch_dynamic<<3) != len )
//		 	++ max_mismatch_dynamic;
//	} else {
//		  max_mismatch_dynamic = ceil( len * kp.mismatch_rate );
//		//max_mismatch_dynamic = floor( len * kp.mismatch_rate );
//	}
//
//	register const char *p = read->seq;
//	for( i=0; i!=len; ++i ) {
//		if( p[pos+i] != kp.adapter_r1[i] ) {
//			if( mis == max_mismatch_dynamic ) return false;
//			++ mis;
//		}
//	}
//
//	return true;
//}

/*
 * use dynamic max_mismatch as the covered size can range from 3 to a large number such as 50,
 * here the maximum mismatch allowed is LEN/4 for read1 + read2
*/
bool check_mismatch_dynamic_PE_C( const CPEREAD *read, unsigned int pos, const ktrim_param &kp ) {
	register unsigned int mis1=0, mis2=0;
	register unsigned int i, len;
	len = read->size - pos;
	if( len > kp.adapter_len )
		len = kp.adapter_len;

	register unsigned int max_mismatch_dynamic;
	// update in v1.1.0: allows the users to set the proportion of mismatches
	// BUT it is highly discouraged
	if( kp.use_default_mismatch ) {
		// each read allows 1/8 mismatches of the total comparable length
		max_mismatch_dynamic = len >> 3;
		// if( (max_mismatch_dynamic<<3) != len )
		// 	++ max_mismatch_dynamic;
	} else {
		max_mismatch_dynamic = floor( len * kp.mismatch_rate );
	}

	// check mismatch for each read
	register const char * p = read->seq1 + pos;
	register const char * q = kp.adapter_r1;
	for( i=0; i!=len; ++i, ++p, ++q ) {
		if( *p != *q ) {
			if( mis1 == max_mismatch_dynamic )
				return false;
			++ mis1;
		}
	}
	p = read->seq2 + pos;
	q = kp.adapter_r2;
	for( i=0; i!=len; ++i, ++p, ++q ) {
		if( *p != *q ) {
			if( mis2 == max_mismatch_dynamic )
				return false;
			++ mis2;
		}
	}

	return true;
}

// update in v1.2: support window check
// p 是要 trim 的序列 total_size是序列处理之前的长度
int get_quality_trim_cycle_se( const char *p, const int total_size, const ktrim_param &kp ) {
	register int i, j, k;
	register int stop = kp.min_length-1;
	for( i=total_size-1; i>=stop; ) {
		if( p[i] >= kp.quality ) {
			k = i - kp.window;
			for( j=i-1; j!=k; --j ) {
				if( j<0 || p[j]<kp.quality ) {
					break;
				}
			}
			if( j == k ) { // find the quality trimming position i
				break;
			} else {	// there is a low-quality base in the middle
				i = j - 1;
			}
		} else {
			-- i;
		}
	}

	if( i >= stop )
		return i + 1; // 返回的是quality-trim之后的序列长度
	else
		return 0; 
}

int get_quality_trim_cycle_pe( const CPEREAD *read, const ktrim_param &kp ) {
	register int i, j, k;
	register int stop = kp.min_length - 1;
	const char *p = read->qual1;
	const char *q = read->qual2;
	for( i=read->size-1; i>=stop; ) {
		if( p[i]>=kp.quality && q[i]>=kp.quality ) {
			k = i - kp.window;
			for( j=i-1; j!=k; --j ) {
				if( j<0 || p[j]<kp.quality || q[j]<kp.quality ) {
					break;
				}
			}
			if( j == k ) { // find the quality trimming position
				break;
			} else {
				i = j - 1;
			}
		} else {
			-- i;
		}
	}
	
	if( i >= stop )
		return i + 1;
	else
		return 0;
}

#endif

