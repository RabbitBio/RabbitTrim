#ifndef ILLUMINA_LONG_CLIPPING_SEQ_H
#define ILLUMINA_LONG_CLIPPING_SEQ_H

#include "IlluminaClippingSeq.h"
#include "../Logger.h"
#include <string>
#include <assert.h>
#include <set>

namespace rabbit
{
    namespace trim
    {
        typedef unsigned long long uint64;
        class IlluminaLongClippingSeq : public IlluminaClippingSeq {
                public:
                    IlluminaLongClippingSeq(rabbit::Logger& logger_, int phred_, std::string seq_, int seedMaxMiss_, int minSequenceLikelihood_, int minSequenceOverlap_);
                    ~IlluminaLongClippingSeq();
                    
                    int packCh(char ch);
                    uint64* packSeqExternal(Reference& rec);
                    uint64 calcSingleMask(int length);
                    float calculateDifferenceQuality(Reference& rec, int overlap, int recOffset);
                    float calculateMaximumRange(float* vals, int valsLen);
                    int readsSeqCompare(Reference& rec);

                    
                    
                private:
                    rabbit::Logger logger;
                    int phred;
                    std::string seq;
                    int seqLen;
                    uint64* pack;
                    uint64* fullPack;
                    int seedMax;
                    int seedMaxMiss;
                    int minSequenceLikelihood;
                    int minSequenceOverlap;


            
        };
    } // namespace trim
} // namespace rabbit

#endif
