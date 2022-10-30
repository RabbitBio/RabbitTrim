#ifndef BASE_COUNT_TRIMMER_H
#define BASE_COUNT_TRIMMER_H

#include "trimmer/Trimmer.h"
#include <string>
#include <climits>
#include <assert.h> 

namespace rabbit
{
    namespace trim
    {
        class BaseCountTrimmer : public Trimmer{
                public:
                    // BaseCountTrimmer(std::string bases, int minCount_ = 0, int maxCount_ = std::INT_MAX){
                    BaseCountTrimmer(std::string bases, int minCount_ = 0, int maxCount_ = (1 << 30)){
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
                    void processRecords(std::vector<Reference>& recs, bool isPair = false, bool isReverse = false);
                private:
                    int minCount;
                    int maxCount;
                    int baseSet;
        };
    } // namespace trim
} // namespace rabbit

#endif
