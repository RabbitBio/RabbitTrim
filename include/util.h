#ifndef _RABBIT_TRIM_UTIL_H
#define _RABBIT_TRIM_UTIL_H

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
#include <stdint.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <vector>
#include <string>
#include <string.h>
#include "FastxChunk.h"
#include "Reference.h"
#include "Formater.h"

namespace rabbit
{
	namespace trim
	{
		namespace util
		{
			bool startsWith(std::string const &value, std::string const &starting);
			bool endsWith(std::string const &value, std::string const &ending);

      std::vector<std::string> split(const std::string& str, const std::string& delim);

      std::string& ClearHeadTailSpace(std::string &str);

			// Determine if the file exists or not
			bool isFileExist(const char* filename);

			// format
			int neoGetLine(rabbit::fq::FastqDataChunk *&chunk, uint64_t &pos, uint64_t &len);

			// int chunkFormat(rabbit::fq::FastqDataChunk *fqDataChunk, CSEREAD *read, bool mHasQuality = true);
			// unsigned int chunkFormat_PE(rabbit::fq::FastqDataPairChunk* chunk,CPEREAD * read, bool mHasQuality = true);


			// 多线程时间统计
			double getTime();

			// extract file names 
			void extractFileNames( const char *str, std::vector<std::string> & Rs );

			//load 1 batch of data, using purely C-style
			// unsigned int load_batch_data_SE_C( FILE * fp, CSEREAD *loadingReads, unsigned int num );
			// unsigned int load_batch_data_SE_GZ( gzFile gfp, CSEREAD *loadingReads, unsigned int num );

			// unsigned int load_batch_data_PE_C( FILE *fq1, FILE *fq2, CPEREAD *loadingReads, unsigned int num );

			// unsigned int load_batch_data_PE_GZ( gzFile gfp1, gzFile gfp2, CPEREAD *loadingReads, unsigned int num );

			/*
			* use dynamic max_mismatch as the covered size can range from 3 to a large number such as 50,
			* here the maximum mismatch allowed is LEN/8
			*/
			// bool check_mismatch_dynamic_SE_C( CSEREAD *read, unsigned int pos, const ktrim_param & kp );
			/*
			* use dynamic max_mismatch as the covered size can range from 3 to a large number such as 50,
			* here the maximum mismatch allowed is LEN/4 for read1 + read2
			*/
			// bool check_mismatch_dynamic_PE_C( const CPEREAD *read, unsigned int pos, const ktrim_param &kp );

			// update in v1.2: support window check
			// p 是要 trim 的序列 total_size是序列处理之前的长度
			// int get_quality_trim_cycle_se( const char *p, const int total_size, const ktrim_param &kp );

			// int get_quality_trim_cycle_pe( const CPEREAD *read, const ktrim_param &kp );
		} // namespace util
		
		
	} // namespace trim
	
	
} // namespace rabbit

#endif

