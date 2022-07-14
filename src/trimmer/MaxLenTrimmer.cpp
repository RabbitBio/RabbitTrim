#include "MaxLenTrimmer.h"

using namespace rabbit;

void MaxLenTrimmer::processOneRecord(Reference& rec){
    int len = rec.length;
    rec.length = len <= maxLen ? len : 0; 
}