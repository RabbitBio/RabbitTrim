#ifndef ILLUMINA_SHORT_CLIPPING_SEQ_H
#define ILLUMINA_SHORT_CLIPPING_SEQ_H
#include "trimmer/IlluminaClippingSeq.h"
#include "Logger.h"
#include <string>
#include <assert.h>
#include <set>
#include <vector>
#include <algorithm>
#include "robin_hood.h"

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
                    float calculateDifferenceQuality(neoReference& rec, int overlap, int recOffset);
                    int readsSeqCompare(neoReference& rec);

                    int packCh(char ch, bool reverse);
                    inline uint64 calcSingleMask(int length);
                    float calculateMaximumRange(float* vals, int valsLen);

                    uint64* packSeqExternal(Reference& rec);
                    float calculateDifferenceQuality(Reference& rec, int overlap, int recOffset);
                    int readsSeqCompare(Reference& rec);
                    unsigned long long cnt = 0;
                    unsigned long long total = 0;
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
