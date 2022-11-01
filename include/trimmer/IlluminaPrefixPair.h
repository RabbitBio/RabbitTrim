#ifndef ILLUMINA_PREFIX_PAIR_H
#define ILLUMINA_PREFIX_PAIR_H

#include <string>
#include "Logger.h"
#include "Reference.h"
#include <assert.h>
namespace rabbit
{
    namespace trim
    {
        typedef unsigned long long uint64;
        class IlluminaPrefixPair{
                public:
                    std::string prefix1;
                    std::string prefix2;
                    
                    IlluminaPrefixPair(std::string prefix1_, std::string prefix2_, rabbit::Logger& logger, int phred_, int minPrefix_, int seedMaxMiss_, int minPalindromeLikelihood_, bool palindromeKeepBoth_);
                    ~IlluminaPrefixPair() = default;
                    
                    
                    uint64 packCh(char ch, bool reverse);
                    uint64* packSeqInternal(std::string seq, bool reverse);
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
                    
        };
    } // namespace trim
} // namespace rabbit

#endif
