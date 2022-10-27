#ifndef SLIDING_WINDOW_TRIMMER_H
#define SLIDING_WINDOW_TRIMMER_H

#include "Trimmer.h"

namespace rabbit
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
            void processRecords(std::vector<Reference&> recs, bool isPair = false, bool isReverse = false);
        private:
            int windowLength;
            float requiredQuality;
            float totalRequiredQuality;
            int phred;
   };
} // namespace rabbit


#endif