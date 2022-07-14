#include "AvgQualTrimmer.h"

using namespace rabbit;

void AvgQualTrimmer::processOneRecord(Reference& rec){
    int len = rec.length;
    int cur_headPos = rec.headPos;
    
    int total = 0;
    for(int i = 0; i < len; i++){
        int qual_val = rec.seq.at(i + cur_headPos) == 'N' ? 0 : rec.quality.at(i + cur_headPos) - phred;
        total += qual_val;
    }

    rec.length = total < qual * len ? 0 : len;

}