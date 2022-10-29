/*
 * common.h
 *
 * This header file records the constants used in Ktrim
 *
 * Author: Kun Sun (sunkun@szbl.ac.cn)
 * Date:   Mar, 2020
 * This program is part of the Ktrim package
 *
 **/

/*
 *  Sequencing model:
 *
 *                *read1 -->
 *  5' adapter - sequence sequence sequence - 3' adapter
 *                                <-- read2*
 *
 *  so read1 may contains 3' adapter, name it adapter_r1,
 *  read2 may contains reversed and ACGT-paired 5' adapter, name it adapter_r2
 *
 **/

#ifndef _KTRIM_COMMON_
#define _KTRIM_COMMON_

#include <string>
#include <vector>
#include <assert.h>
#include <cstdio>
using namespace std;

const char * VERSION = "1.2.2 (Jan 2021)";

/*
 * Change.log
 *
 *
 * 1.2.2 Jan 2021
 *	fix bug when "-o" is NOT present but the program does not quit
 *
 * 1.2.1 Nov 2020
 *	fix bug in multi-file handling
 *
 * 1.2.0 Jun 2020
 *	Multiple performance enhancement including GZ-support
 *	and various optimization in paired-end files
 *
*/

// structure of a READ
const unsigned int MAX_READ_ID    = 128;
const unsigned int MAX_READ_CYCLE = 512;

typedef struct {
	char *id;
	char *seq;
	char *qual; // seq 和 qual 遇到 结束标志数字0 
	unsigned int size; // read 的长度
} CSEREAD;

typedef struct {
	char *id1;
	char *seq1;
	char *qual1;
	unsigned int size;

	char *id2;
	char *seq2;
	char *qual2;
	unsigned int size2;
} CPEREAD;

typedef struct {
	unsigned int *dropped; // trim后丢弃的read的数量
	unsigned int *real_adapter; // 根据seed找到的adapter数量
	unsigned int *tail_adapter; // 根据seed找不到 但尾部符合adapter的read数量
} ktrim_stat;

typedef struct {
	char ** buffer1; // buffer 存储经过quality-trim 和 adapter-trim 后满足条件的read,包括read的id、seq和qual
	char ** buffer2;
	unsigned int *b1stored; // b1stored 表示buffer1上一次写入完成后的位置，初始化是0
	unsigned int *b2stored;
} writeBuffer;


// 设置输出的超大数组
const int CHUNK_NUM = 1 << 7;
const int MEM_PER_CHUNK = 1 << 22; // chunk的大小

typedef struct writeBufferTotal{
	int buffer_id_;
	bool isFinished = false;
	char** buffer1;// buffer 存储经过quality-trim 和 adapter-trim 后满足条件的read,包括read的id、seq和qual
	char** buffer2;
	int chunk_num_ = 0;
	int mem_per_chunk_ = 0;
	struct writeBufferTotal* next1 = NULL;
	struct writeBufferTotal* next2 = NULL;
	writeBufferTotal(){}
	writeBufferTotal(int mem_per_chunk){
		chunk_num_ = CHUNK_NUM; // 使用common.h 中定义的 CHUNK_NUM
		mem_per_chunk_ = mem_per_chunk;
		buffer1 = new char*[chunk_num_];
		buffer2 = NULL;
		for(int i =0; i < chunk_num_; i++){
			buffer1[i] = new char[mem_per_chunk];
			memset(buffer1[i],0,mem_per_chunk);
		}
		next1 = NULL;
		next2 = NULL;
	}
	writeBufferTotal(int chunk_num, int mem_per_chunk){
		chunk_num_ = chunk_num;
		mem_per_chunk_ = mem_per_chunk;
		buffer1 = new char*[chunk_num];
		buffer2 = NULL;
		for(int i =0; i < chunk_num; i++){
			buffer1[i] = new char[mem_per_chunk];
			memset(buffer1[i],0,mem_per_chunk);
		}
		next1 = NULL;
		next2 = NULL;		
	}
	writeBufferTotal(int chunk_num, int mem_per_chunk,int buffer_id,bool isPE){
		assert( isPE == true);
		buffer_id_ = buffer_id;
		chunk_num_ = chunk_num;
		mem_per_chunk_ = mem_per_chunk;
		buffer1 = new char*[chunk_num];
		buffer2 = new char*[chunk_num];
		for(int i =0; i < chunk_num; i++){
			buffer1[i] = new char[mem_per_chunk];
			memset(buffer1[i],0,mem_per_chunk);
			buffer2[i] = new char[mem_per_chunk];
			memset(buffer2[i],0,mem_per_chunk);
		}
		next1 = NULL;
		next2 = NULL;		
	}
	
	void Reset() {
		next1 = NULL;
		next2 = NULL;
		if(buffer2){
			for(int i =0; i < chunk_num_; i++){
				memset(buffer1[i],0,mem_per_chunk_);
				memset(buffer2[i],0,mem_per_chunk_);
			}
		}else{
			for(int i = 0; i < chunk_num_; i++){
				// buffer1[i] = new char[mem_per_chunk_];
				memset(buffer1[i],0,mem_per_chunk_);
			}
		}
	}
	void free(){
		assert(chunk_num_ > 0);
		if(buffer2){
			for (int i = 0; i < chunk_num_; i++){
				delete [] buffer1[i];
				delete [] buffer2[i];
			}
			delete [] buffer1;
			delete [] buffer2;
		}else{
			for (int i = 0; i < chunk_num_; i++){
				delete [] buffer1[i];
			}
			delete [] buffer1;
		}
	}
	// TODO std::atomic_int size; // 如何进行初始化
}writeBufferTotal;

// built-in adapters
const unsigned int MIN_ADAPTER_SIZE = 8;
const unsigned int MAX_ADAPTER_SIZE = 64;
const unsigned int ADAPTER_BUFFER_SIZE = 128;
const unsigned int ADAPTER_INDEX_SIZE = 3;
const unsigned int OFFSET_INDEX3 = 3;
// illumina TruSeq kits adapters
//const char * illumina_adapter_r1 = "AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAG";
//const char * illumina_adapter_r2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT";
//const unsigned int illumina_adapter_len = 33;
const char * illumina_adapter_r1 = "AGATCGGAAGAGC";
const char * illumina_adapter_r2 = illumina_adapter_r1;
const unsigned int illumina_adapter_len = 13;
const char * illumina_index1 = "AGA";			// use first 3 as index
const char * illumina_index2 = illumina_index1;	// use first 3 as index
const char * illumina_index3 = "TCG";			// for single-end data
// Nextera kits and AmpliSeq for Illumina panels
const char * nextera_adapter_r1 = "CTGTCTCTTATACACATCT";
const char * nextera_adapter_r2 = nextera_adapter_r1;
const unsigned int nextera_adapter_len = 19;
const char * nextera_index1 = "CTG";			// use first 3 as index
const char * nextera_index2 = nextera_index1;	// use first 3 as index
const char * nextera_index3 = "TCT";			// for single-end data
// Nextera transposase adapters
const char * transposase_adapter_r1 = "TCGTCGGCAGCGTC";
const char * transposase_adapter_r2 = "GTCTCGTGGGCTCG";
const unsigned int transposase_adapter_len = 14;
const char * transposase_index1 = "TCG";		// use first 3 as index
const char * transposase_index2 = "GTC";		// use first 3 as index
const char * transposase_index3 = "TCG";		// for single-end data
// BGI adapters
//const char * bgi_adapter_r1 = "AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA";
//const char * bgi_adapter_r2 = "AAGTCGGATCGTAGCCATGTCGTTCTGTGAGC";
//const unsigned int bgi_adapter_len = 32;
const char * bgi_adapter_r1 = "AAGTCGGAGGCCAAGCGGTC";
const char * bgi_adapter_r2 = "AAGTCGGATCGTAGCCATGT";
const unsigned int bgi_adapter_len = 20;
const char * bgi_index1 = "AAG";		// use first 3 as index
const char * bgi_index2 = "AAG";		// use first 3 as index
const char * bgi_index3 = "TCG";		// for single-end data

// seed and error configurations
const unsigned int impossible_seed	= 10000;
const unsigned int MAX_READ_LENGTH	= 1024;
const unsigned int MAX_SEED_NUM		= 128;

//configurations for parallelization, which is highly related to memory usage
//but seems to have very minor effect on running time
const int READS_PER_BATCH			= 1 << 16;	// process 256 K reads per batch (for parallelization)
const int BUFFER_SIZE_PER_BATCH_READ= 1 << 26;	// 256 MB buffer for each thread to store FASTQ
const int MEM_SE_READSET = READS_PER_BATCH * (MAX_READ_ID+MAX_READ_CYCLE+MAX_READ_CYCLE);
const int MEM_PE_READSET = READS_PER_BATCH * (MAX_READ_ID+MAX_READ_CYCLE+MAX_READ_CYCLE) * 2;


// enlarge the buffer for single-thread run
//const int READS_PER_BATCH_ST = READS_PER_BATCH << 1;	// process 256 K reads per batch (for parallelization)
//const int BUFFER_SIZE_PER_BATCH_READ_ST = BUFFER_SIZE_PER_BATCH_READ << 1;	// 256 MB buffer for each thread to store FASTQ
const int READS_PER_BATCH_ST			= 1 << 16;	// process 256 K reads for single-thread
const int BUFFER_SIZE_PER_BATCH_READ_ST	= 1 << 26;	// 256 MB buffer for single-thread to store FASTQ

const char FILE_SEPARATOR = ',';

// paramaters related
typedef struct {
	char *FASTQ1, *FASTQ2, *FASTQU;
	char *outpre;

	unsigned int thread;
	unsigned int min_length;
	unsigned int phred;
	unsigned int minqual;
	unsigned int quality; // quality = phred + minqual  即计算出来的read中质量分数的最小值
	unsigned int window;

	const char *seqKit, *seqA, *seqB;
	const char *adapter_r1, *adapter_r2;
	unsigned int adapter_len;
	const char *adapter_index1, *adapter_index2, *adapter_index3;
	// index1 是read1的adapter索引（1-3bp）
	// index2 是read2的adapter索引（1-3bp）
	// index3 是read1的adapter第二个索引（4-6bp）

	bool use_default_mismatch;
	float mismatch_rate;
	
	// trimmomatic params
	char* trim_log;
	bool validatePairing;
	char* trim_steps;
	bool quiet;
	
	RabbitTrimParam(){}
	
} RabbitTrimParam;

const char * param_list = "1:2:U:o:t:k:s:p:q:w:a:b:m:hvl:c:d:e";

// definition of functions
void usage();
void init_param( RabbitTrimParam &rp );
int  process_cmd_param( int argc, char * argv[], RabbitTrimParam &rp );
void print_param( const RabbitTrimParam &rp );
void extractFileNames( const char *str, vector<string> & Rs );

// C-style
unsigned int load_batch_data_PE_C( FILE * fq1, FILE * fq2, CPEREAD *loadingReads, unsigned int num );
bool check_mismatch_dynamic_PE_C( const CPEREAD *read, unsigned int pos, const RabbitTrimParam &rp );
int process_single_thread_PE_C( const RabbitTrimParam &rp );
int process_multi_thread_PE_C(  const RabbitTrimParam &rp );

unsigned int load_batch_data_SE_C( FILE * fp, CSEREAD *loadingReads, unsigned int num );
bool check_mismatch_dynamic_SE_C( const char * p, unsigned int pos, const RabbitTrimParam &rp );
int process_single_thread_SE_C( const RabbitTrimParam &rp );
int process_multi_thread_SE_C(  const RabbitTrimParam &rp );

#endif

