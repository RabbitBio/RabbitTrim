#include "trimmer/MaxLenTrimmer.h"

using namespace rabbit::trim;

void MaxLenTrimmer::processOneRecord(Reference& rec){
    int len = rec.length;
    rec.length = len <= maxLen ? len : 0; 
}
void MaxLenTrimmer::processRecords(std::vector<Reference>& recs, bool isPair, bool isReverse){
    for(Reference& rec : recs){
        processOneRecord(rec);
    }
}

void MaxLenTrimmer::processOneRecord(neoReference& rec){
    int len = rec.lseq;
    rec.lseq = len <= maxLen ? len : 0; 
    rec.lqual = rec.lseq;
}
void MaxLenTrimmer::processRecords(std::vector<neoReference>& recs, bool isPair, bool isReverse){
    for(neoReference& rec : recs){
        processOneRecord(rec);
    }
}
