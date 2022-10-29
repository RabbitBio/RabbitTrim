#include "MinLenTrimmer.h"

using namespace rabbit::trim;

void MinLenTrimmer::processOneRecord(Reference& rec){
    int len = rec.length;
    if(len < minLen) rec.length = 0;
}

void MinLenTrimmer::processRecords(std::vector<Reference>& recs, bool isPair, bool isReverse){
    for(Reference& rec : recs){
        processOneRec(rec);
    }
}