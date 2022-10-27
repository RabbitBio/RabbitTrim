#ifndef TAIL_CROP_TRIMMER_H
#define TAIL_CROP_TRIMMER_H

#include "Trimmer.h"
#include <limits.h>

namespace rabbit{
    class TailCropTrimmer : public Trimmer{
        public: 
            TailCropTrimmer(int bases_, int maxLength_ = std::INT_MAX){
                bases = bases_;
                maxLength = maxLength_;
            }
            ~TailCropTrimmer(){}

            void processOneRecord(Reference& rec);
            void processRecords(std::vector<Reference&> recs, bool isPair = false, bool isReverse = false);
        private:
            int bases;
            int maxLength;

    };
    
}
#endif