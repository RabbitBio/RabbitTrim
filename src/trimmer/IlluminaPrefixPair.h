#ifndef ILLUMINA_PREFIX_PAIR_H
#define ILLUMINA_PREFIX_PAIR_H

#include <string>
#include "../Logger.h"
#include "../io/Reference.h"
#include "../io/Globals.h"
#include <assert.h>
#include <limits.h>
namespace rabbit
{
   class IlluminaPrefixPair{
        public:
            std::string prefix1;
            std::string prefix2;
            
            IlluminaPrefixPair(std::string prefix1_, std::string prefix2_, rabbit::Logger& logger, int phred_, int minPrefix_, int seedMaxMiss_, int minPalindromeLikelihood_, bool palindromeKeepBoth_);
            ~IlluminaPrefixPair() = default;
            
            
            int packCh(char ch, bool reverse);
            uint64* packSeqInternal(std::string seq, bool reverse);
            float calculatePalindromeDifferenceQuality(Reference& rec1, Reference& rec2, int overlap, int skip1, int skip2);
            int palindromeReadsCompare(Reference& rec1, Reference& rec2);
        private:
            constexpr int BASE_A = 1;
            constexpr int BASE_C = 4;
            constexpr int BASE_G = 8;
            constexpr int BASE_T = 2;
            constexpr float LOG10_4 = 0.60206f;
            
            rabbit::Logger logger;
            int phred;
            int prefixLen;
            int seedMaxMiss;
            int seedMax;
            int minPalindromeLikelihood;
            int minPrefix;
            bool palindromeKeepBoth;
            
   } ;
} // namespace rabbit

#endif
