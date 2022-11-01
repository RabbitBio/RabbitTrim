#ifndef MIN_LEN_TRIMMER_H
#define MIN_LEN_TRIMMER_H

#include "trimmer/Trimmer.h"

namespace rabbit
{
    namespace trim
    {
        class MinLenTrimmer : public Trimmer {
            public:
                MinLenTrimmer(int minLen_){
                    minLen = minLen_;
                }
                ~MinLenTrimmer() = default;
            
                void processOneRecord(Reference& rec);
                void processRecords(std::vector<Reference>& recs, bool isPair = false, bool isReverse = false);
                void processOneRecord(neoReference& rec);
                void processRecords(std::vector<neoReference>& recs, bool isPair = false, bool isReverse = false);
            private:
                int minLen;
        };
    } // namespace trim
} // namespace rabbit

#endif
