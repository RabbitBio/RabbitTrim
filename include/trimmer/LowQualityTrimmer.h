#ifndef LOW_QUALITY_TRIMMER_H
#define LOW_QUALITY_TRIMMER_H
#include "trimmer/Trimmer.h"

namespace rabbit
{
    namespace trim
    {
        class LowQualityTrimmer : public Trimmer
        {
        private:
            int quality;
            int window;
            int phred;

        public:
            LowQualityTrimmer(int quality_, int window_, int phred_);
            ~LowQualityTrimmer();
            void processOneRecord(Reference& rec);
            void processRecords(std::vector<Reference>& recs, bool isPair, bool isReverse);
            void processOneRecord(neoReference& rec);
            void processRecords(std::vector<neoReference>& recs, bool isPair, bool isReverse);
        };
        
        LowQualityTrimmer::LowQualityTrimmer(int quality_, int window_, int phred_)
        {
            quality = quality_;
            window = window_;
            phred = phred_;
        }
        
        LowQualityTrimmer::~LowQualityTrimmer()
        {
            
        }
        
    } // namespace trim
    
    
} // namespace rabbit

#endif
