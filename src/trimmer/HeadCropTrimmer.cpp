#include "trimmer/HeadCropTrimmer.h"

using namespace rabbit::trim;

void HeadCropTrimmer::processOneRecord(Reference& rec){
    int len = rec.length;
    int toTrim;
    toTrim = len - bases > maxLength ? len - maxLength : bases;
    toTrim = len < toTrim ? len : toTrim;
    rec.headPos += toTrim;
    rec.length = len - toTrim;
    
}

void HeadCropTrimmer::processRecords(std::vector<Reference>& recs, bool isPair, bool isReverse){
    for(Reference& rec : recs){
        processOneRecord(rec);
    }
}
        

void HeadCropTrimmer::processOneRecord(neoReference& rec){
    int len = rec.lseq;
    int toTrim;
    toTrim = len - bases > maxLength ? len - maxLength : bases;
    toTrim = len < toTrim ? len : toTrim;
    rec.pseq += toTrim;
    rec.pqual += toTrim;
    rec.lseq = len - toTrim;
    rec.lqual = rec.lseq;
    
}

void HeadCropTrimmer::processRecords(std::vector<neoReference>& recs, bool isPair, bool isReverse){
    for(neoReference& rec : recs){
        processOneRecord(rec);
    }
}
void HeadCropTrimmer::processRecords(std::vector<neoReference>& recs, int threadId, bool isPair, bool isReverse){
    for(neoReference& rec : recs){
        processOneRecord(rec);
    }
}
