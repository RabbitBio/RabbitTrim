#ifndef ILLUMINA_CLIPPING_TRIMMER_H
#define ILLUMINA_CLIPPING_TRIMMER_H

#include "trimmer/Trimmer.h"
#include "Logger.h"
#include "trimmer/IlluminaPrefixPair.h"
#include "trimmer/IlluminaClippingSeq.h"
#include "trimmer/IlluminaShortClippingSeq.h"
#include "trimmer/IlluminaMediumClippingSeq.h"
#include "trimmer/IlluminaLongClippingSeq.h"
#include <fstream>
#include <vector>
#include <map>
#include <assert.h>
#include "util.h"
#include <limits.h>
namespace rabbit
{
    namespace trim
    {
        class IlluminaClippingTrimmer : public Trimmer{
                public:
                    IlluminaClippingTrimmer(rabbit::Logger& logger_, int phred_, std::string fastaAdapterFile, int seedMaxMiss_, int minPalindromeLikelihood_, int minSequenceLikelihood_, int minPrefix_ = 8, bool palindromeKeepBoth_ = false);
                    ~IlluminaClippingTrimmer();
                    
                    void processOneRecord(Reference& rec);
                    void processSingleRecord(Reference& rec, bool isReverse = false);
                    void processPairRecord(Reference& rec1, Reference& rec2);
                    void processRecords(std::vector<Reference>& recs, bool isPair = false, bool isReverse = false);

                    void processOneRecord(neoReference& rec);
                    void processSingleRecord(neoReference& rec, bool isReverse = false);
                    void processPairRecord(neoReference& rec1, neoReference& rec2);
                    void processRecords(std::vector<neoReference>& recs, bool isPair = false, bool isReverse = false);
                private:
                    const std::string PREFIX = "Prefix";
                    const std::string SUFFIX_F = "/1";
                    const std::string SUFFIX_R = "/2";
                    const int INTERLEAVE = 4;
                    const float LOG10_4 = 0.60206f;
                    const int BASE_A = 1;
                    const int BASE_C = 4;
                    const int BASE_G = 8;
                    const int BASE_T = 2;

                    rabbit::Logger logger;
                    int phred;
                    int seedMaxMiss;
                    int minPalindromeLikelihood;
                    int minSequenceLikelihood;
                    int minSequenceOverlap;
                    int minPrefix;
                    bool palindromeKeepBoth;
                    
                    std::vector<IlluminaPrefixPair*> prefixPairs;
                    std::vector<IlluminaClippingSeq*> forwardSeqs;
                    std::vector<IlluminaClippingSeq*> reverseSeqs;
                    std::vector<IlluminaClippingSeq*> commonSeqs;
        };
    } // namespace trim
} // namespace rabbit

#endif
