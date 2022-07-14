#include "MinLenTrimmer.h"

using namespace rabbit;

void MinLenTrimmer::processOneRecord(Reference& rec){
    int len = rec.length;
    if(len < minLen) rec.length = 0;
}