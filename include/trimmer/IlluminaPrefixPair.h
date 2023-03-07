#ifndef ILLUMINA_PREFIX_PAIR_H
#define ILLUMINA_PREFIX_PAIR_H

#include <string>
#include "Logger.h"
#include "Reference.h"
#include <assert.h>
#include "param.h"
#include <algorithm>
namespace rabbit
{
    namespace trim
    {
        typedef unsigned long long uint64;
        class IlluminaPrefixPair{
                public:
                    std::string prefix1;
                    std::string prefix2;
                    
                    IlluminaPrefixPair(std::string prefix1_, std::string prefix2_, rabbit::Logger& logger, int phred_, int minPrefix_, int seedMaxMiss_, int minPalindromeLikelihood_, bool palindromeKeepBoth_, int consumerNum_ = 1);
                    ~IlluminaPrefixPair();
                    
                    
                    uint64 packCh(char ch, bool reverse);
                    uint64* packSeqInternal(std::string seq, bool reverse);

                    
                    uint64* packSeqInternalForward(neoReference& rec, int threadId);
                    uint64* packSeqInternalReverse(neoReference& rec, int threadId);
                    int palindromeReadsCompare(neoReference& rec1, neoReference& rec2, int threadId);
                    float calculatePalindromeDifferenceQuality(neoReference& rec1, neoReference& rec2, int overlap, int skip1, int skip2, int threadId);


                    uint64* packSeqInternal(neoReference& rec, bool reverse);
                    float calculatePalindromeDifferenceQuality(neoReference& rec1, neoReference& rec2, int overlap, int skip1, int skip2);
                    int palindromeReadsCompare(neoReference& rec1, neoReference& rec2);

                    uint64* packSeqInternal(Reference& rec, bool reverse);
                    float calculatePalindromeDifferenceQuality(Reference& rec1, Reference& rec2, int overlap, int skip1, int skip2);
                    int palindromeReadsCompare(Reference& rec1, Reference& rec2);
                private:
                    const uint64 BASE_A = 1;
                    const uint64 BASE_C = 4;
                    const uint64 BASE_G = 8;
                    const uint64 BASE_T = 2;
                    const float LOG10_4 = 0.60206f;
                    
                    rabbit::Logger logger;
                    int phred;
                    int prefixLen;
                    int seedMaxMiss;
                    int seedMax;
                    int minPalindromeLikelihood;
                    int minPrefix;
                    bool palindromeKeepBoth;
                    int consumerNum;
                    uint64* forwardPacks;
                    uint64* reversePacks;
                    float* likelihoodArr;
                    
        };
    } // namespace trim
} // namespace rabbit

#endif
