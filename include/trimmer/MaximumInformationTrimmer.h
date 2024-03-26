#ifndef MAXIMUM_INFORMATION_TRIMMER_H
#define MAXIMUM_INFORMATION_TRIMMER_H

#include "trimmer/Trimmer.h"
#include "param.h"
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

                MaximumInformationTrimmer(int parLength_, float strictness_, int phred_, int consumerNum_);
                ~MaximumInformationTrimmer();

                void processOneRecord(Reference& rec);
                void processRecords(std::vector<Reference>& recs, bool isPair = false, bool isReverse = false);
                void processOneRecord(neoReference& rec);
                void processOneRecord(neoReference& rec, int threadId);
                void processRecords(std::vector<neoReference>& recs, bool isPair = false, bool isReverse = false);
                void processRecords(std::vector<neoReference>& recs, int threadId, bool isPair = false, bool isReverse = false);

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
                
                char* quals;
#if defined __SSE2__ && defined __SSE4_1__ && defined TRIM_USE_VEC
                char* all_N;
                char* phred_arr;
                char* max_qual;
#endif

        };
    } // namespace trim
    
} // namespace rabbit


#endif
