#ifndef HEAD_CROP_TRIMMER_H
#define HEAD_CROP_TRIMMER_H

#include "Trimmer.h"
#include <limits.h>
namespace rabbit
{
    namespace trim
    {
        class HeadCropTrimmer : public Trimmer{
                public:
                    // HeadCropTrimmer(int bases_, int maxLength_ = std::INT_MAX){
                    HeadCropTrimmer(int bases_, int maxLength_ = (1 << 30)){
                        bases = bases_;
                        maxLength = maxLength_;    
                    }
                    ~HeadCropTrimmer(){}
                    
                    void processOneRecord(Reference& rec);
                    void processRecords(std::vector<Reference>& recs, bool isPair = false, bool isReverse = false);

                private:
                    int bases;
                    int maxLength; 
        };
    } // namespace trim
} // namespace rabbit

#endif
