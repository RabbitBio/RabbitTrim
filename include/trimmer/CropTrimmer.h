#ifndef CROP_TRIMMER_H
#define CROP_TRIMMER_H

#include "trimmer/Trimmer.h"
namespace rabbit
{
    namespace trim
    {
        class CropTrimmer : public Trimmer{
            public:
                CropTrimmer(int len_){
                    len = len_;
                }
                
                ~CropTrimmer(){}
                
                void processOneRecord(Reference& rec);
                void processRecords(std::vector<Reference>& recs, bool isPair = false, bool isReverse = false);
                void processOneRecord(neoReference& rec);
                void processRecords(std::vector<neoReference>& recs, bool isPair = false, bool isReverse = false);
                
            private:
                int len;
        };
    } // namespace trim
} // namespace rabbit

#endif
