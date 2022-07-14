#include "HeadCropTrimmer.h"

using namespace rabbit;

void HeadCropTrimmer::processOneRecord(Reference& rec){
    int len = rec.length;
    int toTrim;
    toTrim = len - bases > maxLength ? len - maxLength : bases;
    toTrim = len < toTrim ? len : toTrim;
    rec.headPos += toTrim;
    rec.length = len - toTrim;
    
}