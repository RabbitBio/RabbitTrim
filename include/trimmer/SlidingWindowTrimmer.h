#ifndef SLIDING_WINDOW_TRIMMER_H
#define SLIDING_WINDOW_TRIMMER_H

#include "trimmer/Trimmer.h"

namespace rabbit
{
    namespace trim
    {
        class SlidingWindowTrimmer : public Trimmer {
                public:
                    SlidingWindowTrimmer(int windowLength_, float requiredQuality_, int phred_){
                        windowLength = windowLength_;
                        requiredQuality = requiredQuality_;
                        totalRequiredQuality = requiredQuality * windowLength;
                        phred = phred_;
                    }
                    ~SlidingWindowTrimmer(){}
                    
                    void processOneRecord(Reference& rec);
                    void processRecords(std::vector<Reference>& recs, bool isPair = false, bool isReverse = false);
                    void processOneRecord(neoReference& rec);
                    void processRecords(std::vector<neoReference>& recs, bool isPair = false, bool isReverse = false);
                    void processRecords(std::vector<neoReference>& recs, int threadId, bool isPair = false, bool isReverse = false);
                private:
                    int windowLength;
                    float requiredQuality;
                    float totalRequiredQuality;
                    int phred;
        };
    } // namespace trim
} // namespace rabbit


#endif
