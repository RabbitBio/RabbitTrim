#ifndef CROP_TRIMMER_H
#define CROP_TRIMMER_H

#include "Trimmer.h"
namespace rabbit
{
    class CropTrimmer : public Trimmer{
        public:
            CropTrimmer(int len_){
                len = len_;
            }
            
            ~CropTrimmer(){}
            
            void processOneRecord(Reference& rec);
            
        private:
            int len;
    };
    
} // namespace rabbit

#endif
