#ifndef AVG_QUAL_TRIMMER_H
#define AVG_QUAL_TRIMMER_H

#include "Trimmer.h"

namespace rabbit
{
   class AvgQualTrimmer : public Trimmer{
        public:
            AvgQualTrimmer(int qual_, int phred_){
                qual = qual_;
                phred = phred_;
            }
            ~AvgQualTrimmer() = default;

            void processOneRecord(Reference& rec);
        private:
            int qual;
            int phred;
   } ;
} // namespace rabbit

#endif