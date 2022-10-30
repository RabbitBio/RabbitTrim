#ifndef MAXIMUM_INFORMATION_TRIMMER_H
#define MAXIMUM_INFORMATION_TRIMMER_H

#include "trimmer/Trimmer.h"
#include <cmath>
#include <limits>

namespace rabbit
{
    namespace trim
    {
        typedef long long int64;
        class MaximumInformationTrimmer : public Trimmer{
            public:
                const int LONGEST_READ = 1000;
                const int MAXQUAL = 60;

                MaximumInformationTrimmer(int parLength_, float strictness_);
                ~MaximumInformationTrimmer() = default;

                void processOneRecord(Reference& rec);
                void processRecords(std::vector<Reference>& recs, bool isPair = false, bool isReverse = false);

            private:
                int phred;
                int parLength; // target length
                float strictness;
                
                double* lengthScoreTmp;
                double* qualProbTmp;
                int64* lengthScore;
                int64* qualProb;
                
                double calcNormalization(double* arr, int arrLength, int margin);
                int64* normalize(double* arr, int arrLength, double ratio);

        };
    } // namespace trim
    
} // namespace rabbit


#endif
