#ifndef MAX_LEN_TRIMMER_H
#define MAX_LEN_TRIMMER_H

#include "trimmer/Trimmer.h"

namespace rabbit
{
    namespace trim
    {
        class MaxLenTrimmer : public Trimmer{
                public:
                    MaxLenTrimmer(int maxLen_){
                        maxLen = maxLen_;
                    }
                    ~MaxLenTrimmer() = default;

                    void processOneRecord(Reference& rec);
                    void processRecords(std::vector<Reference>& recs, bool isPair = false, bool isReverse = false);
                    void processOneRecord(neoReference& rec);
                    void processRecords(std::vector<neoReference>& recs, bool isPair = false, bool isReverse = false);
                    void processRecords(std::vector<neoReference>& recs, int threadId, bool isPair = false, bool isReverse = false);
                private: 
                    int maxLen;
        };
    } // namespace trim
} // namespace rabbit


#endif
