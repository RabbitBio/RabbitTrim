#ifndef ILLUMINA_CLIPPING_SEQ_H
#define ILLUMINA_CLIPPING_SEQ_H
#include <string>
#include "../io/Globals.h"
#include "../io/Reference.h"

namespace rabbit
{
   class IlluminaClippingSeq {
        public:
            virtual int readsSeqCompare(Reference& rec) = 0;
            virtual float calculateDifferenceQuality(Reference& rec, int overlap, int recOffset) = 0;
            // virtual float calculateMaximumRange(float* vals, int valsLen) = 0;
            


            ~IlluminaClippingSeq() = default;
        private:
            constexpr int BASE_A = 1;
            constexpr int BASE_C = 4;
            constexpr int BASE_G = 8;
            constexpr int BASE_T = 2;
            constexpr uint64 MASK_VAL = 0xffffffffffffffffULL;
            constexpr float LOG10_4 = 0.60206f;
            constexpr int INTERLEAVE = 4;
            
   } ;
} // namespace rabbit

#endif