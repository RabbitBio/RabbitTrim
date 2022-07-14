#ifndef TO_PHREAD_64_TRIMMER_H
#define TO_PHREAD_64_TRIMMER_H

#include "Trimmer.h"

namespace rabbit
{
   class ToPhred64Trimmer : public Trimmer {
        public:
            ToPhred64Trimmer(int phred_){
                phred = phred_;
            }
            ~ToPhred64Trimmer() = default;
            
            void processOneRecord(Reference& rec);
        private:
            int phred;
   };
} // namespace rabbit


#endif