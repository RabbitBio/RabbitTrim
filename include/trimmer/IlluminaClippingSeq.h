#ifndef ILLUMINA_CLIPPING_SEQ_H
#define ILLUMINA_CLIPPING_SEQ_H
#include <string>
#include "Globals.h"
#include "Reference.h"

namespace rabbit
{
    namespace trim
    {
        typedef unsigned long long uint64;
        class IlluminaClippingSeq {
                public:
                    virtual int readsSeqCompare(Reference& rec) = 0;
                    virtual float calculateDifferenceQuality(Reference& rec, int overlap, int recOffset) = 0;
                    // virtual float calculateMaximumRange(float* vals, int valsLen) = 0;
                    ~IlluminaClippingSeq() = default;
                protected:
                    const int BASE_A = 1;
                    const int BASE_C = 4;
                    const int BASE_G = 8;
                    const int BASE_T = 2;
                    const uint64 MASK_VAL = 0xffffffffffffffffULL;
                    const float LOG10_4 = 0.60206f;
                    const int INTERLEAVE = 4;
                    
        };
    } // namespace trim
} // namespace rabbit

#endif
