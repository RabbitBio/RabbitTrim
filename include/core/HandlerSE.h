#ifndef HANDLER_SE_H
#define HANDLER_SE_H
#include <vector>
#include <functional>
#include <mutex>
#include <atomic>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <thread>
#include <ctime>
#include <stdlib.h>
#include <memory.h>
#include <omp.h>
#include "common.h"
#include "PairingValidator.h"
#include "Logger.h"
#include "param.h"
#include "FastxChunk.h"
#include "trimlog/TrimStat.h"
#include "trimmer/TrimmerFactory.h"

namespace rabbit
{
    namespace trim
    {
        typedef rabbit::core::TDataQueue<rabbit::fq::FastqDataChunk> FastqDataChunkQueue;
        int  process_se(rabbit::trim::RabbitTrimParam& rp, rabbit::Logger &logger);
        int  producer_se_task(rabbit::trim::RabbitTrimParam& rp, rabbit::Logger& logger, rabbit::fq::FastqDataPool* fastqPool, FastqDataChunkQueue& dq);
        void consumer_se_task(rabbit::trim::RabbitTrimParam& rp, rabbit::fq::FastqDataPool *fastqPool, FastqDataChunkQueue &dq, rabbit::trim::TrimStat& rstats, std::vector<rabbit::trim::Trimmer*>& trimmers);
    } // namespace trim
    
} // namespace rabbit

#endif

