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
#include "param.h"

namespace rabbit
{
    namespace trim
    {
        // seq.length < 16
        class IlluminaShortClippingSeq : public IlluminaClippingSeq {
                public:
                    IlluminaShortClippingSeq(rabbit::Logger& logger_, int phred_, std::string seq_, int seedMaxMiss_, int minSequenceLikelihood_, int minSequenceOverlap_, int consumerNum_);
                    ~IlluminaShortClippingSeq();
                    
                    uint64* packSeqExternal(neoReference& rec);
                    float calculateDifferenceQuality(neoReference& rec, int overlap, int recOffset);
                    int readsSeqCompare(neoReference& rec);

                    uint64* packSeqExternal(neoReference& rec, int threadId);
                    float calculateDifferenceQuality(neoReference& rec, int overlap, int recOffset, int threadId);
                    int readsSeqCompare(neoReference& rec, int threadId);

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
                    int consumerNum;
                    uint64 curSize;
                    void* worker_buffer;
                    uint64* recPacks; // storage the rec pack for each thread
                    float* likelihood;
                    void reAllocateBuffer(uint64 recLen);

#if defined __SSE2__ && defined __AVX__ && defined __AVX2__
                    char* seq_str; // store adapter and extra 15 'B'
                    float* awards; // all elements are LOG10_4
                    char* phred_arr;
                    char* all_N;
                    float* divide_arr;
#endif

        };
    } // namespace trim
} // namespace rabbit

#endif
