#ifndef MIN_LEN_TRIMMER_H
#define MIN_LEN_TRIMMER_H

#include "Trimmer.h"

namespace rabbit
{
    class MinLenTrimmer : public Trimmer {
        public:
            MinLenTrimmer(int minLen_){
                minLen = minLen_;
            }
            ~MinLenTrimmer() = default;
        
            void processOneRecord(Reference& rec);
        private:
            int minLen;
    };
} // namespace rabbit

#endif
