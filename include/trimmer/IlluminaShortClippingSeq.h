#ifndef ILLUMINA_SHORT_CLIPPING_SEQ_H
#define ILLUMINA_SHORT_CLIPPING_SEQ_H
#include "trimmer/IlluminaClippingSeq.h"
#include "Logger.h"
#include <string>
#include <assert.h>
#include <set>
#include <vector>

namespace rabbit
{
    namespace trim
    {
        // seq.length < 16
        class IlluminaShortClippingSeq : public IlluminaClippingSeq {
                public:
                    IlluminaShortClippingSeq(rabbit::Logger& logger_, int phred_, std::string seq_, int seedMaxMiss_, int minSequenceLikelihood_, int minSequenceOverlap_);
                    ~IlluminaShortClippingSeq();
                    
                    uint64* packSeqExternal(neoReference& rec);
                    uint64* packSeqExternal(Reference& rec);
                    int packCh(char ch, bool reverse);
                    uint64 calcSingleMask(int length);
                    float calculateMaximumRange(float* vals, int valsLen);

                    float calculateDifferenceQuality(Reference& rec, int overlap, int recOffset);
                    float calculateDifferenceQuality(Reference& rec, int overlap, int recOffset);
                    int readsSeqCompare(Reference& rec);
                private:
                    rabbit::Logger& logger;
                    int phred;
                    std::string seq; // fasta seq
                    int seqLen;
                    uint64* pack;
                    uint64 mask;
                    int seedMax;
                    int seedMaxMiss;
                    int minSequenceLikelihood;
                    int minSequenceOverlap;

        };
    } // namespace trim
} // namespace rabbit

#endif
