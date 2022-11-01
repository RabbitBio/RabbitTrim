#ifndef ILLUMINA_LONG_CLIPPING_SEQ_H
#define ILLUMINA_LONG_CLIPPING_SEQ_H

#include "trimmer/IlluminaClippingSeq.h"
#include "Logger.h"
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
                    uint64 calcSingleMask(int length);
                    float calculateMaximumRange(float* vals, int valsLen);

                    uint64* packSeqExternal(Reference& rec);
                    float calculateDifferenceQuality(Reference& rec, int overlap, int recOffset);
                    int readsSeqCompare(Reference& rec);

                    uint64* packSeqExternal(neoReference& rec);
                    float calculateDifferenceQuality(neoReference& rec, int overlap, int recOffset);
                    int readsSeqCompare(neoReference& rec);

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
