#ifndef TO_PHRED_33_TRIMMER_H
#define TO_PHRED_33_TRIMMER_H

#include "Trimmer.h"
namespace rabbit
{
   class ToPhred33Trimmer : public Trimmer {
        public:
            ToPhred33Trimmer() = default;
            ToPhred33Trimmer(int phred_){
                phred = phred_;
            }
            ~ToPhred33Trimmer() = default;
            
            void processOneRecord(Reference& rec);
            void processRecords(std::vector<Reference&> recs, bool isPair = false, bool isReverse = false);
        
        private:
            int phred;
   };
} // namespace rabbit

#endif