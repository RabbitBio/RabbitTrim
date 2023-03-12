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

namespace rabbit
{
    namespace trim
    {
        typedef rabbit::core::TDataQueue<rabbit::fq::FastqDataChunk> FastqDataChunkQueue;
        typedef rabbit::core::TDataQueue<rabbit::trim::WriterBuffer> WriterDataQueue;
        typedef rabbit::core::TDataQueue<rabbit::log::TrimLogBuffer> TrimLogDataQueue;

        // output data pool
        typedef rabbit::core::TDataPool<rabbit::trim::WriterBuffer> WriterBufferDataPool;

        int  process_se(rabbit::trim::RabbitTrimParam& rp, rabbit::Logger &logger);
        int  producer_se_task(rabbit::trim::RabbitTrimParam& rp, rabbit::Logger& logger, rabbit::fq::FastqDataPool* fastqPool, FastqDataChunkQueue& dq);
        void consumer_se_task(rabbit::trim::RabbitTrimParam& rp, rabbit::fq::FastqDataPool *fastqPool, FastqDataChunkQueue &dq, WriterBufferDataPool* wbDataPool, WriterDataQueue& dq2, TrimLogDataQueue& dq3, rabbit::log::TrimStat& rstats, const std::vector<rabbit::trim::Trimmer*>& trimmers, int threadId, std::atomic_ullong& atomic_next_id);
        void writer_se_task(rabbit::trim::RabbitTrimParam& rp, WriterBufferDataPool* wbDataPool, WriterDataQueue& dq2, rabbit::Logger& logger); 
        void trimlog_se_task(rabbit::trim::RabbitTrimParam& rp, TrimLogDataQueue& dq, rabbit::Logger& logger);

        // pigz
        void pigzer_se_task(rabbit::trim::RabbitTrimParam& rp, WriterBufferDataPool* wbDataPool, WriterDataQueue& dq2, std::pair<char*, int>& pigzLast);
    } // namespace trim

} // namespace rabbit

#endif

