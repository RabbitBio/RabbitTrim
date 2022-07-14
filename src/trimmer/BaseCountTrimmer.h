#ifndef BASE_COUNT_TRIMMER_H
#define BASE_COUNT_TRIMMER_H

#include "Trimmer.h"
#include <string>
#include <limits.h>
#include <assert.h> 

namespace rabbit
{
   class BaseCountTrimmer : public Trimmer{
        public:
            BaseCountTrimmer(std::string bases, int minCount_ = 0; int maxCount_ = std::INT_MAX){
                minCount = minCount_;
                maxCount = maxCount_;
                
                int nums = bases.size();
                baseSet = 0;
                for(int i = 0; i < nums; i++){
                    int pos = bases.at(i) - 'A';
                    assert(pos >= 0);
                    baseSet |= 1 << pos ; 
                }

            }
            ~BaseCountTrimmer() = default;

            void processOneRecord(Reference& rec);
        private:
            int minCount;
            int maxCount;
            int baseSet;
   } ;
} // namespace rabbit

#endif