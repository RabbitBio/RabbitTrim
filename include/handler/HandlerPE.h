#ifndef HANDLER_PE_H
#define HANDLER_PE_H
#include <vector>
#include <functional>
#include <mutex>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <thread>
#include <ctime>
#include <stdlib.h>
#include "PairingValidator.h"
#include "Logger.h"
#include "param.h"
#include "util.h"
#include "FastxChunk.h"
#include "trimlog/TrimStat.h"
#include "trimmer/TrimmerFactory.h"
#include "Reference.h"
#include "handler/WriterBuffer.h"

namespace rabbit
{
    namespace trim
    {
        typedef rabbit::core::TDataQueue<rabbit::fq::FastqDataPairChunk> FastqDataPairChunkQueue;
        int  process_pe(rabbit::trim::RabbitTrimParam& rp, rabbit::Logger &logger);
        int  producer_pe_task(rabbit::trim::RabbitTrimParam& rp, rabbit::Logger& logger, rabbit::fq::FastqDataPool* fastqPool, FastqDataPairChunkQueue& dq);
        void consumer_pe_task(rabbit::trim::RabbitTrimParam& rp, rabbit::fq::FastqDataPool *fastqPool, FastqDataPairChunkQueue &dq, PEWriterDataQueue& dq2, rabbit::trim::TrimStat& rstats, std::vector<rabbit::trim::Trimmer*>& trimmers);
        void writer_pe_task(rabbit::trim::RabbitTrimParam& rp, PEWriterDataQueue& dq2, rabbit::Logger& logger);
    } // namespace trim
    
} // namespace rabbit

#endif
