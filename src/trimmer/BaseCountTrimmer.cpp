#include "trimmer/BaseCountTrimmer.h"

using namespace rabbit::trim;

void BaseCountTrimmer::processOneRecord(Reference& rec){
    int count = 0;
    int len = rec.length;
    int cur_headPos = rec.headPos;
    for(int i = 0; i < len; i++){
        int pos = rec.seq.at(i + cur_headPos) - 'A';
        assert(pos >= 0) ;
        if((1 << pos) & baseSet) count++; 
    }
    if(count < minCount || count > maxCount) rec.length = 0;
}

void BaseCountTrimmer::processRecords(std::vector<Reference>& recs, bool isPair, bool isReverse){
    for(Reference& rec : recs){
        processOneRecord(rec);
    }
}

void BaseCountTrimmer::processOneRecord(neoReference& rec){
    int count = 0;
    int len = rec.lseq;
    char* seq = (char*)(rec.base + rec.pseq);
    for(int i = 0; i < len; i++){
        int pos = seq[i] - 'A';
        assert(pos >= 0);
        if((1 << pos) & baseSet) count++; 
    }
    if(count < minCount || count > maxCount) rec.lseq = 0;
    rec.lqual = rec.lseq;
}

void BaseCountTrimmer::processRecords(std::vector<neoReference>& recs, bool isPair, bool isReverse){
    for(neoReference& rec : recs){
        processOneRecord(rec);
    }
}
void BaseCountTrimmer::processRecords(std::vector<neoReference>& recs, int threadId, bool isPair, bool isReverse){
    for(neoReference& rec : recs){
        processOneRecord(rec);
    }
}
