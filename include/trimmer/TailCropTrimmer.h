#ifndef TAIL_CROP_TRIMMER_H
#define TAIL_CROP_TRIMMER_H

#include "trimmer/Trimmer.h"
#include <limits.h>

namespace rabbit{
    namespace trim
    {
        class TailCropTrimmer : public Trimmer{
            public: 
                // TailCropTrimmer(int bases_, int maxLength_ = std::INT_MAX){
                TailCropTrimmer(int bases_, int maxLength_ = 1 << 30){
                    bases = bases_;
                    maxLength = maxLength_;
                }
                ~TailCropTrimmer(){}

                void processOneRecord(Reference& rec);
                void processRecords(std::vector<Reference>& recs, bool isPair = false, bool isReverse = false);
                void processOneRecord(neoReference& rec);
                void processRecords(std::vector<neoReference>& recs, bool isPair = false, bool isReverse = false);
                void processRecords(std::vector<neoReference>& recs, int threadId, bool isPair = false, bool isReverse = false);
            private:
                int bases;
                int maxLength;

        };
    } // namespace trim
}
#endif
