#ifndef HANDLER_PE_H
#define HANDLER_PE_H
#include <vector>
#include <string>
#include <functional>
#include <mutex>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <thread>
#include <ctime>
#include <stdlib.h>
#include <atomic>
#include "PairingValidator.h"
#include "Logger.h"
#include "param.h"
#include "util.h"
#include "FastxChunk.h"
#include "trimlog/TrimStat.h"
#include "trimlog/TrimLogRecord.h"
#include "trimmer/TrimmerFactory.h"
#include "Reference.h"
#include "handler/WriterBuffer.h"
#include "handler/HandlerSE.h"
#include "pragzip.h"


namespace rabbit
{
    namespace trim
    {
        typedef rabbit::core::TDataQueue<rabbit::fq::FastqDataPairChunk> FastqDataPairChunkQueue;
        typedef rabbit::core::TDataQueue<rabbit::trim::WriterBuffer> WriterBufferDataQueue;
        // typedef rabbit::core::TDataQueue<rabbit::log::TrimLogBuffer> TrimLogDataQueue;
        // typedef rabbit::core::TDataQueue<rabbit::trim::PEWriterBuffer> PEWriterDataQueue;
        
        // output data pool
        typedef rabbit::core::TDataPool<rabbit::trim::WriterBuffer> WriterBufferDataPool;

        int process_pe(rabbit::trim::RabbitTrimParam& rp, rabbit::Logger &logger);
        int producer_pe_task(rabbit::trim::RabbitTrimParam& rp, rabbit::Logger& logger, rabbit::fq::FastqDataPool* fastqPool, PragzipQueue& pragzipQueue1, PragzipQueue& pragzipQueue2, FastqDataPairChunkQueue& dq, FastqDataPairChunkQueue& phredQueue, std::atomic_bool& isDeterminedPhred);
        void consumer_pe_task(rabbit::trim::RabbitTrimParam& rp, rabbit::fq::FastqDataPool *fastqPool, FastqDataPairChunkQueue &dq, WriterBufferDataPool* wbDataPool, WriterBufferDataQueue& wbQueue1, WriterBufferDataQueue& wbQueue2,  WriterBufferDataQueue& wbQueue3, WriterBufferDataQueue& wbQueue4, WriterBufferDataQueue& logQueue, rabbit::log::TrimStat& rstats, int threadId, std::atomic_ullong& atomic_next_id, rabbit::Logger& logger);
        void writer_pe_task(rabbit::trim::RabbitTrimParam& rp, WriterBufferDataPool* wbDataPool, WriterBufferDataQueue& wbQueue, std::string out_file, rabbit::Logger& logger);

        void consumer_pe_task2(rabbit::trim::RabbitTrimParam& rp, rabbit::fq::FastqDataPool *fastqPool, FastqDataPairChunkQueue &dq, WriterBufferDataPool* wbDataPool, WriterBufferDataQueue& wbQueue1, WriterBufferDataQueue& wbQueue2, rabbit::log::TrimStat& rstats, std::atomic_ullong& atomic_next_id, rabbit::Logger& logger);

        // pigzer
        void pigzer_pe_task(rabbit::trim::RabbitTrimParam& rp, WriterBufferDataPool* wbDataPool, WriterBufferDataQueue& wbQueue, std::pair<char*, int>& pigzLast, std::string out_file); // out_file is the name without ".gz"

        // pragzip 
        void pragzip_pe_task(rabbit::trim::RabbitTrimParam& rp, PragzipQueue& pragzipQueue, std::string in_file);
    } // namespace trim
    
} // namespace rabbit

#endif
