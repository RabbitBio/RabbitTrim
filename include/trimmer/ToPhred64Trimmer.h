#ifndef TO_PHREAD_64_TRIMMER_H
#define TO_PHREAD_64_TRIMMER_H

#include "trimmer/Trimmer.h"
#include <cassert>

namespace rabbit
{
    namespace trim
    {
        class ToPhred64Trimmer : public Trimmer {
                public:
                    ToPhred64Trimmer(int phred_){
                        phred = phred_;
                    }
                    ~ToPhred64Trimmer() = default;
                    
                    void processOneRecord(Reference& rec);
                    void processRecords(std::vector<Reference>& recs, bool isPair = false, bool isReverse = false);
                    void processOneRecord(neoReference& rec);
                    void processRecords(std::vector<neoReference>& recs, bool isPair = false, bool isReverse = false);
                private:
                    int phred;
        };
    } // namespace trim
} // namespace rabbit


#endif
