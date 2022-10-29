#include "CropTrimmer.h"

using namespace rabbit::trim;

void CropTrimmer::processOneRecord(Reference& rec){
    int cur_len = rec.length;
    cur_len = cur_len > len ? len : cur_len;
    rec.length = cur_len;
}

void CropTrimmer::processRecords(std::vector<Reference>& recs, bool isPair, bool isReverse){
    for(Reference& rec : recs){
        processOneRec(rec);
    }
}