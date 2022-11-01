#ifndef TRAILING_TRIMMER_H
#define TRAILING_TRIMMER_H

#include "trimmer/Trimmer.h"
#include <string>

namespace rabbit
{
    namespace trim
    {
        class TrailingTrimmer : public Trimmer{
            public:
                TrailingTrimmer(int qual_, int phred_){
                    qual = qual_;
                    phred = phred_;
                }
                ~TrailingTrimmer(){}
                
                void processOneRecord(Reference& rec);
                void processRecords(std::vector<Reference>& recs, bool isPair = false, bool isReverse = false);
                void processOneRecord(neoReference& rec);
                void processRecords(std::vector<neoReference>& recs, bool isPair = false, bool isReverse = false);
            private:
                int qual;
                int phred;
        };
    } // namespace trim
    
} // namespace rabbit


#endif
