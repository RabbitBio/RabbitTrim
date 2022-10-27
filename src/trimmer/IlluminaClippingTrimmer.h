#ifndef ILLUMINA_CLIPPING_TRIMMER_H
#define ILLUMINA_CLIPPING_TRIMMER_H

#include "Trimmer.h"
#include "../Logger.h"
#include "IlluminaPrefixPair.h"
#include "IlluminaClippingSeq.h"
#include "IlluminaShortClippingSeq.h"
#include "IlluminaMediumClippingSeq.h"
#include "IlluminaLongClippingSeq.h"
#include <fstream>
#include <vector>
#include <map>
#include <assert.h>
#include "../util.h"
#include <limits.h>
namespace rabbit
{
   class IlluminaClippingTrimmer : public Trimmer{
        public:
            IlluminaClippingTrimmer(rabbit::Logger& logger_, int phred_, std::string fastaAdapterFile, int seedMaxMiss_, int minPalindromeLikelihood_, int minSequenceLikelihood_, int minPrefix_ = 8, int palindromeKeepBoth_ = false);
            ~IlluminaClippingTrimmer();
            
            void processOneRecord(Reference& rec);
            void processSingleRecord(Reference& rec, bool isReverse = false);
            void processPairRecord(Reference& rec1, Reference& rec2);
            void processRecords(std::vector<Reference&> recs, bool isPair = false, bool isReverse = false);
        private:
            constexpr std::string PREFIX = "Prefix";
            constexpr std::string SUFFIX_F = "/1";
            constexpr std::string SUFFIX_R = "/2";
            constexpr int INTERLEAVE = 4;
            constexpr float LOG10_4 = 0.60206f;
            constexpr int BASE_A = 1;
            constexpr int BASE_C = 4;
            constexpr int BASE_G = 8;
            constexpr int BASE_T = 2;

            rabbit::Logger logger;
            int phred;
            int seedMaxMiss;
            int minPalindromeLikelihood;
            int minSequenceLikelihood;
            int minSequenceOverlap;
            int minPrefix;
            bool palindromeKeepBoth;
            
            std::vector<IlluminaPrefixPair&> prefixPairs;
            std::vector<IlluminaClippingSeq&> forwardSeqs;
            std::vector<IlluminaClippingSeq&> reverseSeqs;
            std::vector<IlluminaClippingSeq&> commonSeqs;
            

            
            
   } ;
} // namespace rabbit

#endif
