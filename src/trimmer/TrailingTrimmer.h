#ifndef TRAILING_TRIMMER_H
#define TRAILING_TRIMMER_H

#include "Trimmer.h"

namespace rabbit
{
    class TrailingTrimmer : public Trimmer{
        public:
            TrailingTrimmer(int qual_, int phred_){
                qual = qual_;
                phred = phred_;
            }
            ~TrailingTrimmer(){}
            
            void processOneRecord(Reference& rec);
            void processRecords(std::vector<Reference&> recs, bool isPair = false, bool isReverse = false);
        private:
            int qual;
            int phred;
    };
    
} // namespace rabbit


#endif
