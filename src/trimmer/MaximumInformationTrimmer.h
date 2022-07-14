#ifndef MAXIMUM_INFORMATION_TRIMMER_H
#define MAXIMUM_INFORMATION_TRIMMER_H

#include "Trimmer.h"
#include <cmath>

namespace rabbit
{
    class MaximumInformationTrimmer : public Trimmer{
        public:
            constexpr int LONGEST_READ = 1000;
            constexpr int MAXQUAL = 60;

            MaximumInformationTrimmer(int parLength_, float strictness_);
            ~MaximumInformationTrimmer() = default;

            void processOneRecord(Reference& rec);

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
    
} // namespace rabbit


#endif