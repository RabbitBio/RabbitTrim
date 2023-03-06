#ifndef TO_PHRED_33_TRIMMER_H
#define TO_PHRED_33_TRIMMER_H

#include "trimmer/Trimmer.h"
#include <cassert>
namespace rabbit
{
    namespace trim
    {
        class ToPhred33Trimmer : public Trimmer {
                public:
                    ToPhred33Trimmer() = default;
                    ToPhred33Trimmer(int phred_){
                        phred = phred_;
                    }
                    ~ToPhred33Trimmer() = default;
                    
                    void processOneRecord(Reference& rec);
                    void processRecords(std::vector<Reference>& recs, bool isPair = false, bool isReverse = false);
                    void processOneRecord(neoReference& rec);
                    void processRecords(std::vector<neoReference>& recs, bool isPair = false, bool isReverse = false);
                    void processRecords(std::vector<neoReference>& recs, int threadId, bool isPair = false, bool isReverse = false);
                
                private:
                    int phred;
        };
    } // namespace trim
} // namespace rabbit

#endif
