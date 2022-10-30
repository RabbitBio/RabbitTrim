#include "trimmer/HeadCropTrimmer.h"

using namespace rabbit::trim;

void HeadCropTrimmer::processOneRecord(Reference& rec){
    int len = rec.length;
    int toTrim;
    toTrim = len - bases > maxLength ? len - maxLength : bases;
    toTrim = len < toTrim ? len : toTrim;
    rec.headPos += toTrim;
    rec.length = len - toTrim;
    
}

void HeadCropTrimmer::processRecords(std::vector<Reference>& recs, bool isPair, bool isReverse){
    for(Reference& rec : recs){
        processOneRecord(rec);
    }
}
        
