#include "LeadingTrimmer.h"
#include <string>

using namespace rabbit;
    
void LeadingTrimmer::processOneRecord(Reference& rec){
    std::string quality = rec.quality;
    std::string seq = rec.seq;
    int len = rec.length;
    int cur_headPos = rec.headPos;
    for(int i = 0; i < len; i++){
        // N 的质量分数当做0
        char ch = seq.at(i + cur_headPos);
        int qual_val = quality.at(i + cur_headPos) - phred;
        if(ch != 'N' && qual_val >= qual){
            // rec.seq = seq.substr(i);
            // rec.quality = quality.substr(i);
            rec.headPos += i;
            rec.length = len - i;
            return;
        }
        
    }
    // The entire sequence quality score is not met
    // rec.quality = "";
    // rec.seq = "";
    rec.headPos += len;
    rec.length = 0;
}

void LeadingTrimmer::processRecords(std::vector<Reference&> recs, bool isPair, bool isReverse){
    for(Reference& rec : recs){
        processOneRec(rec);
    }
}