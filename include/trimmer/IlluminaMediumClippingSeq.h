#ifndef ILLUMINA_MEDIUM_CLIPPING_SEQ
#define ILLUMINA_MEDIUM_CLIPPING_SEQ
#include "trimmer/IlluminaClippingSeq.h"
#include "Logger.h"
#include <string>
#include <set>
#include <limits.h>
#include <cassert>

namespace rabbit
{
    namespace trim
    {
        // 16 <= seq.length() < 24
        class IlluminaMediumClippingSeq : public IlluminaClippingSeq{
                public:
                    IlluminaMediumClippingSeq(rabbit::Logger& logger_, int phred_, std::string seq_, int seedMaxMiss_, int minSequenceLikelihood_, int minSequenceOverlap_);
                    ~IlluminaMediumClippingSeq();

                    uint64* packSeqExternal(Reference& rec);
                    float calculateDifferenceQuality(Reference& rec, int overlap, int recOffset);
                    int readsSeqCompare(Reference& rec);

                    uint64* packSeqExternal(neoReference& rec);
                    float calculateDifferenceQuality(neoReference& rec, int overlap, int recOffset);
                    int readsSeqCompare(neoReference& rec);

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
                    
        };
    } // namespace trim
} // namespace rabbit

#endif
