#ifndef AVG_QUAL_TRIMMER_H
#define AVG_QUAL_TRIMMER_H

#include "trimmer/Trimmer.h"

namespace rabbit
{
    namespace trim
    {
        class AvgQualTrimmer : public Trimmer{
                public:
                    AvgQualTrimmer(int qual_, int phred_){
                        qual = qual_;
                        phred = phred_;
                    }
                    ~AvgQualTrimmer() = default;

                    void processOneRecord(Reference& rec);
                    void processRecords(std::vector<Reference>& recs, bool isPair = false, bool isReverse = false);
                private:
                    int qual;
                    int phred;
        };
    } // namespace trim
} // namespace rabbit

#endif