#include "trimmer/CropTrimmer.h"

using namespace rabbit::trim;

void CropTrimmer::processOneRecord(Reference& rec){
    int cur_len = rec.length;
    cur_len = cur_len > len ? len : cur_len;
    rec.length = cur_len;
}

void CropTrimmer::processRecords(std::vector<Reference>& recs, bool isPair, bool isReverse){
    for(Reference& rec : recs){
        processOneRecord(rec);
    }
}

void CropTrimmer::processOneRecord(neoReference& rec){
    int cur_len = rec.lseq;
    cur_len = cur_len > len ? len : cur_len;
    rec.lseq = cur_len;
    rec.lqual = rec.lqual;
}

void CropTrimmer::processRecords(std::vector<neoReference>& recs, bool isPair, bool isReverse){
    for(neoReference& rec : recs){
        processOneRecord(rec);
    }
}
void CropTrimmer::processRecords(std::vector<neoReference>& recs, int threadId, bool isPair, bool isReverse){
    for(neoReference& rec : recs){
        processOneRecord(rec);
    }
}
