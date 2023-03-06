#include "trimmer/MinLenTrimmer.h"

using namespace rabbit::trim;

void MinLenTrimmer::processOneRecord(Reference& rec){
    int len = rec.length;
    if(len < minLen) rec.length = 0;
}

void MinLenTrimmer::processRecords(std::vector<Reference>& recs, bool isPair, bool isReverse){
    for(Reference& rec : recs){
        processOneRecord(rec);
    }
}
void MinLenTrimmer::processOneRecord(neoReference& rec){
    int len = rec.lseq;
    if(len < minLen) rec.lseq = 0;
    rec.lqual = rec.lseq;
}

void MinLenTrimmer::processRecords(std::vector<neoReference>& recs, bool isPair, bool isReverse){
    for(neoReference& rec : recs){
        processOneRecord(rec);
    }
}
void MinLenTrimmer::processRecords(std::vector<neoReference>& recs, int threadId, bool isPair, bool isReverse){
    for(neoReference& rec : recs){
        processOneRecord(rec);
    }
}
