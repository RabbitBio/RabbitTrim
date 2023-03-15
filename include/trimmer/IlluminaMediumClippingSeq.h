#ifndef ILLUMINA_MEDIUM_CLIPPING_SEQ
#define ILLUMINA_MEDIUM_CLIPPING_SEQ
#include "trimmer/IlluminaClippingSeq.h"
#include "Logger.h"
#include <string>
#include <set>
#include <limits.h>
#include <cassert>
#include "param.h"

namespace rabbit
{
    namespace trim
    {
        // 16 <= seq.length() < 24
        class IlluminaMediumClippingSeq : public IlluminaClippingSeq{
                public:
                    IlluminaMediumClippingSeq(rabbit::Logger& logger_, int phred_, std::string seq_, int seedMaxMiss_, int minSequenceLikelihood_, int minSequenceOverlap_, int consumerNum_);
                    ~IlluminaMediumClippingSeq();

                    uint64* packSeqExternal(Reference& rec);
                    float calculateDifferenceQuality(Reference& rec, int overlap, int recOffset);
                    int readsSeqCompare(Reference& rec);

                    uint64* packSeqExternal(neoReference& rec);
                    float calculateDifferenceQuality(neoReference& rec, int overlap, int recOffset);
                    int readsSeqCompare(neoReference& rec);

                    uint64* packSeqExternal(neoReference& rec, int threadId);
                    float calculateDifferenceQuality(neoReference& rec, int overlap, int recOffset, int threadId);
                    int readsSeqCompare(neoReference& rec, int threadId);

                    int packCh(char ch);
                    uint64 calcSingleMask(int length);
                    float calculateMaximumRange(float* vals, int valsLen);

                private:
                    rabbit::Logger logger;
                    int phred;
                    std::string seq;
                    int seqLen;
                    uint64* pack;
                    int seedMax;
                    int seedMaxMiss;
                    int minSequenceLikelihood;
                    int minSequenceOverlap;
                    int consumerNum;
                    uint64* recPacks;
                    float* likelihoodTotal;
                    
#if defined __SSE2__ && defined __AVX__ && defined __AVX2__ && defined TRIM_USE_VEC
                    char* seq_str;
                    float* awards;
                    char* phred_arr;
                    char* all_N;
                    float* divide_arr;
#endif
                    
                    
        };
    } // namespace trim
} // namespace rabbit

#endif
