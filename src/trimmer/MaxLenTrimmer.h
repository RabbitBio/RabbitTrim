#ifndef MAX_LEN_TRIMMER_H
#define MAX_LEN_TRIMMER_H

#include "Trimmer.h"

namespace rabbit
{
   class MaxLenTrimmer : public Trimmer{
        public:
            MaxLenTrimmer(int maxLen_){
                maxLen = maxLen_;
            }
            ~MaxLenTrimmer() = default;

            void processOneRecord(Reference& rec);
        private: 
            int maxLen;
   } ;
} // namespace rabbit


#endif