#ifndef HANDLER_SE_H
#define HANDLER_SE_H
#include <vector>
#include <functional>
#include <mutex>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <thread>
#include <ctime>
#include <stdlib.h>
#include <string.h>
#include <atomic>
#include "util.h"
#include "PairingValidator.h"
#include "Logger.h"
#include "param.h"
#include "FastxChunk.h"
#include "trimlog/TrimStat.h"
#include "trimlog/TrimLogRecord.h"
#include "trimmer/TrimmerFactory.h"
#include "handler/WriterBuffer.h"
#include "pragzip.h"

namespace rabbit
{
    namespace trim
    {
        typedef rabbit::core::TDataQueue<rabbit::fq::FastqDataChunk> FastqDataChunkQueue;
        typedef rabbit::core::TDataQueue<rabbit::trim::WriterBuffer> WriterBufferDataQueue;
        // typedef rabbit::core::TDataQueue<rabbit::log::TrimLogBuffer> TrimLogDataQueue;

        // output data pool
        typedef rabbit::core::TDataPool<rabbit::trim::WriterBuffer> WriterBufferDataPool;

        int  process_se(rabbit::trim::RabbitTrimParam& rp, rabbit::Logger &logger);
        int  producer_se_task(rabbit::trim::RabbitTrimParam& rp, rabbit::Logger& logger, rabbit::fq::FastqDataPool* fastqPool,  PragzipQueue& pragzipQueue, FastqDataChunkQueue& dq, FastqDataChunkQueue& phredQueue, std::atomic_bool& isDeterminedPhred);
        void consumer_se_task(rabbit::trim::RabbitTrimParam& rp, rabbit::fq::FastqDataPool *fastqPool, FastqDataChunkQueue &dq, WriterBufferDataPool* wbDataPool, WriterBufferDataQueue& wbQueue, WriterBufferDataQueue& logQueue, rabbit::log::TrimStat& rstats, const std::vector<rabbit::trim::Trimmer*>& trimmers, int threadId, std::atomic_ullong& atomic_next_id);
        void writer_se_task(rabbit::trim::RabbitTrimParam& rp, WriterBufferDataPool* wbDataPool, WriterBufferDataQueue& wbQueue, std::string out_file, rabbit::Logger& logger); 

        // pigz
        void pigzer_se_task(rabbit::trim::RabbitTrimParam& rp, WriterBufferDataPool* wbDataPool, WriterBufferDataQueue& wbQueue, std::pair<char*, int>& pigzLast, std::string out_file); // out_file is the output file name without ".gz"
        // pragzip
        void pragzip_se_task(rabbit::trim::RabbitTrimParam& rp, PragzipQueue& pragzip); 
    } // namespace trim

} // namespace rabbit

#endif

