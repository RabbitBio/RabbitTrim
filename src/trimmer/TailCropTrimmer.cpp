#include "TailCropTrimmer.h"

using namespace rabbit::trim;

void TailCropTrimmer::processOneRecord(Reference& rec){
    int len = rec.length;
    int toTrim;
    toTrim = len - bases > maxLength ? len - maxLength : bases;
    toTrim = len < toTrim ? len : toTrim;
    rec.length = len - toTrim;
}

void TailCropTrimmer::processRecords(std::vector<Reference>& recs, bool isPair, bool isReverse){
    for(Reference& rec : recs){
        processOneRec(rec);
    }
}