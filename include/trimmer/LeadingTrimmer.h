#ifndef LEADING_TRIMMER_H
#define LEADING_TRIMMER_H
#include "trimmer/Trimmer.h"
#include "Reference.h"

namespace rabbit
{
    namespace trim
    {
        class LeadingTrimmer : public Trimmer{
            private:
                int qual;    
                int phred;
            public:
                LeadingTrimmer(int qual_, int phred_){
                    qual = qual_;
                    phred = phred_;
                }
                ~LeadingTrimmer(){}

                // Trimmer Interface
                void processOneRecord(Reference& rec);
                void processRecords(std::vector<Reference>& recs, bool isPair = false, bool isReverse = false);
                void processOneRecord(neoReference& rec);
                void processRecords(std::vector<neoReference>& recs, bool isPair = false, bool isReverse = false);
                void processRecords(std::vector<neoReference>& recs, int threadId, bool isPair = false, bool isReverse = false);
                
        };
    } // namespace trim
    
    
} // namespace rabbit

#endif
