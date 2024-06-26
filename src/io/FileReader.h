#include "Buffer.h"
#include "FastxChunk.h"
#include "utils.h"
#include <iostream>
#include <string>
#include <cstring>
#include "Reference.h"
#include <assert.h>

#if defined(USE_IGZIP)
#include "igzip_lib.h"
#else
#include <zlib.h>  //support gziped files, functional but inefficient
#endif

#if defined(_WIN32)
#define _CRT_SECURE_NO_WARNINGS
#pragma warning(disable : 4996)  // D_SCL_SECURE
#pragma warning(disable : 4244)  // conversion uint64 to uint32
//# pragma warning(disable : 4267)
#define FOPEN fopen
#define FDOPEN fdopen
#define FSEEK _fseeki64
#define FTELL _ftelli64
#define FCLOSE fclose
#elif __APPLE__  // Apple by default suport 64 bit file operations (Darwin 10.5+)
#define FOPEN fopen
#define FDOPEN fdopen
#define FSEEK fseek
#define FTELL ftell
#define FCLOSE fclose
#else
#if !defined(_LARGEFILE_SOURCE)
#define _LARGEFILE_SOURCE
#if !defined(_LARGEFILE64_SOURCE)
#define _LARGEFILE64_SOURCE
#endif
#endif
#if defined(_FILE_OFFSET_BITS) && (_FILE_OFFSET_BITS != 64)
#undef _FILE_OFFSET_BITS
#endif
#if !defined(_FILE_OFFSET_BITS)
#define _FILE_OFFSET_BITS 64
#endif
#define FOPEN fopen64
#define FDOPEN fdopen
#define FSEEK fseeko64
#define FTELL ftello64
#define FCLOSE fclose
#endif

namespace rabbit{
	class FileReader{
	private:
		static const uint32 IGZIP_IN_BUF_SIZE = 1 << 22; // 4M gziped file onece fetch
		static const uint32 GZIP_HEADER_BYTES_REQ = 1<<16;
	public:
    typedef rabbit::core::TDataQueue<std::pair<char*, int> > PragzipQueue;
		FileReader(const std::string &fileName_, bool isZipped){
			if(ends_with(fileName_, ".gz") || isZipped) {
#if defined(USE_IGZIP)
				// printf("using igzip!!\n");
				mFile = fopen(fileName_.c_str(), "rb");
				if (mFile == NULL){
					throw RioException(
						("Can not open file to read: " + fileName_).c_str());  
				}
				mIgInbuf = new unsigned char[IGZIP_IN_BUF_SIZE];
				isal_gzip_header_init(&mIgzipHeader);
				isal_inflate_init(&mStream);
				mStream.crc_flag = ISAL_GZIP_NO_HDR_VER;

				mStream.next_in = mIgInbuf;
				mStream.avail_in = fread(mStream.next_in, 1, IGZIP_IN_BUF_SIZE, mFile);

				int ret =0;
				ret = isal_read_gzip_header(&mStream, &mIgzipHeader);
				if( ret != ISAL_DECOMP_OK ){
					cerr << "error invalid gzip header found: " << fileName_ << endl;
					if(mFile != NULL){
						fclose(mFile);
					}
					exit(-1);
				}
#else
				mZipFile = gzopen(fileName_.c_str(), "r");
				if (mZipFile == NULL) {
				  throw RioException(("Can not open file to read: "+fileName_).c_str());  
				}
				gzrewind(mZipFile);
#endif			
				this->isZipped = true;
			}else {
				mFile = FOPEN(fileName_.c_str(), "rb");
				if (fileName_ != "") {
					mFile = FOPEN(fileName_.c_str(), "rb");
					if (mFile == NULL)
						throw RioException("Can not open file to read: ");  
				}
				if (mFile == NULL) {
					throw RioException(
						("Can not open file to read: " + fileName_).c_str());  
				}
			}
		}

		FileReader(int fd, bool isZipped = false){
			if (isZipped) {
#if defined(USE_IGZIP)
				mFile = fdopen(fd, "r");
				if (mFile == NULL)
				{
					cerr << "Can not open file id to read: " + to_string(fd) << endl;
				}
				mIgInbuf = new unsigned char[IGZIP_IN_BUF_SIZE];
				isal_gzip_header_init(&mIgzipHeader);
				isal_inflate_init(&mStream);
				mStream.crc_flag = ISAL_GZIP_NO_HDR_VER;

				mStream.next_in = mIgInbuf;
				mStream.avail_in = fread(mStream.next_in, 1, IGZIP_IN_BUF_SIZE, mFile);

				int ret = 0;
				ret = isal_read_gzip_header(&mStream, &mIgzipHeader);
				if (ret != ISAL_DECOMP_OK)
				{
					cerr << "error! invalid gzip header found!" << endl;
					if (mFile != NULL)
					{
						fclose(mFile);
					}
					exit(-1);
				}
#else
				mZipFile = gzdopen(fd, "r");
				gzrewind(mZipFile);
#endif
				this->isZipped = true;
			} else {
				mFile = FDOPEN(fd, "rb");
				if (fd != -1) {
					mFile = FDOPEN(fd, "rb");
					if (mFile == NULL)
						throw RioException("Can not open file to read: ");  
				}
				if (mFile == NULL) {
					throw RioException("Can not open file to read: ");  
				}
			}
		}

#if	defined(USE_IGZIP)
		int64 igzip_read(FILE* zipFile, byte *memory_, size_t size_){
			uint64_t offset = 0;
			int ret = 0;
			do{
				if(mStream.avail_in == 0){
					mStream.next_in = mIgInbuf;
					mStream.avail_in = fread(mStream.next_in, 1, IGZIP_IN_BUF_SIZE, zipFile);
				}
				mStream.next_out = memory_ + offset;
				mStream.avail_out = size_ - offset;
				if(isal_inflate(&mStream) != ISAL_DECOMP_OK){
					cerr << "decompress error" << endl;
					return -1;
				}
				offset = (mStream.next_out - memory_);

				if(mStream.block_state == ISAL_BLOCK_FINISH) {
					if(feof(mFile) && mStream.avail_in == 0){
						return offset;
					}
					// a new block begins
					if (mStream.avail_in == 0) {
						isal_inflate_reset(&mStream);
						mStream.next_in = mIgInbuf;
						mStream.avail_in = fread(mStream.next_in, 1, IGZIP_IN_BUF_SIZE, mFile);
						//mGzipInputUsedBytes += mStream.avail_in;
					} else if (mStream.avail_in >= GZIP_HEADER_BYTES_REQ){
						unsigned char* old_next_in = mStream.next_in;
						size_t old_avail_in = mStream.avail_in;
						isal_inflate_reset(&mStream);
						mStream.avail_in = old_avail_in;
						mStream.next_in = old_next_in;
					} else {
						size_t old_avail_in = mStream.avail_in;
						memmove(mIgInbuf, mStream.next_in, mStream.avail_in);
						size_t readed = 0;
						if(!feof(mFile)){
							readed = fread(mIgInbuf + mStream.avail_in, 1, IGZIP_IN_BUF_SIZE - mStream.avail_in, mFile);
						}
						isal_inflate_reset(&mStream);
						mStream.next_in = mIgInbuf;
						mStream.avail_in = old_avail_in + readed;;
					}
					ret = isal_read_gzip_header(&mStream, &mIgzipHeader);
					if (ret != ISAL_DECOMP_OK) {
						cerr << "igzip: invalid gzip header found: " << mStream.avail_in <<  " : " << mStream.avail_out << endl;
						exit(-1);
					}
				}
			}while(mStream.avail_out > 0);
			assert(offset <= size_);
			return offset;
		}
#endif

		int64 Read(byte *memory_, uint64 size_) {
			if (isZipped) {
#if defined(USE_IGZIP)			
				int64 n = igzip_read(mFile, memory_, size_);
#else
				int64 n = gzread(mZipFile, memory_, size_);
#endif
				if (n == -1) std::cerr << "Error to read gzip file" << std::endl;
				return n;
			} else {
				int64 n = fread(memory_, 1, size_, mFile);
				return n;
			}
		}
    
    // Add for pragzip
        int64 Read(byte *buf, uint64 len, PragzipQueue& pragzipQueue, std::pair<char *, int> &L,int tag) {
          static int64 chunk_id;
          static int now_L_size[2];
          now_L_size[0] = 1 << 20;
          now_L_size[1] = 1 << 20;
          std::pair<char *, int>* now;
          int64 ret;
          int64 got = 0;
          //printf("now producer get data from pugz queue\n");
          if (L.second > 0) {
            if (L.second >= len) {
              memcpy(buf, L.first, len);
              char *tmp = new char[L.second - len];
              memcpy(tmp, L.first + len, L.second - len);
              memcpy(L.first, tmp, L.second - len);
              delete[] tmp;
              L.second = L.second - len;
              ret = len;
              buf += ret;
              len -= ret;
              got += ret;
              return got;
            } else {
              memcpy(buf, L.first, L.second);
              ret = L.second;
              L.second = 0;
              buf += ret;
              len -= ret;
              got += ret;
            }
          }
          // bool overWhile = false;
          while (len > 0) {
            // while (Q->try_dequeue(now) == 0) {
            //   if (Q->size_approx() == 0 && *done == 1) {
            //     ret = 0;
            //     overWhile = true;
            //     break;
            //   }
            //   usleep(1000);
            //   //printf("producer waiting pugz...\n");
            // }
            if(pragzipQueue.Pop(chunk_id, now) == false)
            {
              ret = 0;
              break;
            }
            // if (overWhile) {
            //   ret = 0;
            //   break;
            // }
            if(now->second > now_L_size[tag]) {
              char* tmp_L = new char[now->second];
              memcpy(tmp_L, L.first, L.second);
              delete[] L.first;
              now_L_size[tag] = now->second;
              L.first = tmp_L;
            }
            if (now->second <= len) {
              memcpy(buf, now->first, now->second);
              delete[] now->first;
              ret = now->second;
            } else {
              int move_last = now->second - len;
              memcpy(buf, now->first, len);
              memcpy(L.first, now->first + len, move_last);
              L.second = move_last;
              delete[] now->first;
              ret = len;
            }

            buf += ret;
            len -= ret;
            got += ret;
          }

          return got;
        }
        /// True means no need to call Read() function
		bool FinishRead(){
			if(isZipped){
#if defined(USE_IGZIP)			
				return feof(mFile) && (mStream.avail_in == 0);
#else
				return gzeof(mZipFile);
#endif			
			}else{
				return feof(mFile);
			}
		}
		bool Eof() const {
			return eof;
			//if(eof) return eof;
			//return feof(mFile);
		}
		void setEof(){
			eof = true;
		}

		~FileReader(){
			if(mIgInbuf != NULL) delete mIgInbuf;
			if(mFile != NULL){
				fclose(mFile);
				mFile = NULL;
			}
			if(mZipFile != NULL){
				gzclose(mZipFile);
			}
		}

	private:
		FILE* mFile = NULL;
		gzFile mZipFile = NULL;
		//igzip usage
		unsigned char *mIgInbuf = NULL;
#if defined(USE_IGZIP)	
		isal_gzip_header mIgzipHeader;
		inflate_state mStream;
#endif	
		bool isZipped = false;
		bool eof = false;
	};

} //namespace rabbit
