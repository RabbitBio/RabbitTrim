#include "trimmer/TailCropTrimmer.h"

using namespace rabbit::trim;

void TailCropTrimmer::processOneRecord(Reference& rec){
    int len = rec.length;
    int toTrim;
    toTrim = len - bases > maxLength ? len - maxLength : bases;
    toTrim = len < toTrim ? len : toTrim;
    rec.length = len - toTrim;
}

void TailCropTrimmer::processRecords(std::vector<Reference>& recs, bool isPair, bool isReverse){
    for(Reference& rec : recs){
        processOneRecord(rec);
    }
}

void TailCropTrimmer::processOneRecord(neoReference& rec){
    int len = rec.lseq;
    int toTrim;
    toTrim = len - bases > maxLength ? len - maxLength : bases;
    toTrim = len < toTrim ? len : toTrim;
    rec.lseq = len - toTrim;
    rec.lqual = rec.lseq; 
}

void TailCropTrimmer::processRecords(std::vector<neoReference>& recs, bool isPair, bool isReverse){
    for(neoReference& rec : recs){
        processOneRecord(rec);
    }
}
void TailCropTrimmer::processRecords(std::vector<neoReference>& recs, int threadId, bool isPair, bool isReverse){
    for(neoReference& rec : recs){
        processOneRecord(rec);
    }
}
