#ifndef ILLUMINA_CLIPPING_TRIMMER_H
#define ILLUMINA_CLIPPING_TRIMMER_H

#include "Trimmer.h"
#include "../Logger.h"
#include <fstream>
namespace rabbit
{
   class IlluminaClippingTrimmer : public Trimmer{
        public:
            IlluminaClippingTrimmer(rabbit::Logger& logger_, std::string fastaAdapterFile, int seedMaxMiss_, int minPalindromeLikelihood_, int minSequenceLikelihood_, int minPrefix_ = 1, int palindromeKeepBoth_ = false);
            ~IlluminaClippingTrimmer();
            
            void processOneRecord(Reference& rec);
        private:
            constexpr std::string PREFIX = "Prefix";
            constexpr std::string SUFFIX_F = "/1";
            constexpr std::string SUFFIX_R = "/2";
            constexpr int INTERLEAVE = 4;
            constexpr float LOG10_4 = 0.60206f;

            rabbit::Logger logger;
            int seedMaxMiss;
            int minPalindromeLikelihood;
            int minSequenceLikelihood;
            int minSequenceOverlap;
            int minPrefix;
            bool palindromeKeepBoth;
            
            

            
            
   } ;
} // namespace rabbit

#endif
