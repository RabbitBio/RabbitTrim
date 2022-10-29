#ifndef HANDLER_PE_H
#define HANDLER_PE_H
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
#include "io/FastxChunk.h"
#include "trimlog/TrimStat.h"
#include "trimmer/TrimmerFactory.h"

namespace rabbit
{
    namespace trim
    {
        typedef rabbit::core::TDataQueue<rabbit::fq::FastqDataPairChunk> FastqDataPairChunkQueue;
        int  process_pe(rabbit::RabbitTrimParam& rp, rabbit::Logger &logger);
        int  producer_pe_task(rabbit::RabbitTrimParam& rp, rabbit::Logger& logger, rabbit::fq::FastqDataPool* fastqPool, FastqDataPairChunkQueue& dq);
        void consumer_pe_task(rabbit::RabbitTrimParam& rp, rabbit::fq::FastqDataPool *fastqPool, FastqDataPairChunkQueue &dq, rabbit::TrimStat& rstats, std::vector<rabbit::trim::Trimmer*>& trimmers);
    } // namespace trim
    
} // namespace rabbit

#endif
