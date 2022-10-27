#include "ToPhred33Trimmer.h"

using namespace rabbit;

void ToPhred33Trimmer::processOneRecord(Reference& rec){
    if(phred == 33) return;
    int len = rec.length;
    int cur_headPos = rec.headPos;
    char* new_quality = new char[len + 1];
    for(int i = 0; i < len; i++){
        char qual_ch = (char)(rec.quality.at(i + cur_headPos) - 31);
        new_quality[i] = qual_ch; 
    }
    new_quality[len] = '\0';
    std::string qualityTmp = std::string(new_quality);
    rec.quality = qualityTmp; 
}

void ToPhred33Trimmer::processRecords(std::vector<Reference&> recs, bool isPair, bool isReverse){
    for(Reference& rec : recs){
        processOneRec(rec);
    }
}