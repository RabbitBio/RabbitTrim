#include "trimmer/AvgQualTrimmer.h"

using namespace rabbit::trim;

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


void AvgQualTrimmer::processRecords(std::vector<Reference>& recs, bool isPair, bool isReverse){
    for(Reference& rec : recs){
        processOneRecord(rec);
    }
}


void AvgQualTrimmer::processOneRecord(neoReference& rec){
    int len = rec.lseq;
    char* rec_seq = (char*)(rec.base +  rec.pseq);
    char* rec_qual = (char*)(rec.base + rec.pqual);
    
    int total = 0;
    for(int i = 0; i < len; i++){
        int qual_val = (rec_seq[i] == 'N' ? 0 : rec_qual[i] - phred);
        total += qual_val;
    }
    rec.lseq = ((total < qual * len) ? 0 : len);
    rec.lqual = rec.lseq;

}


void AvgQualTrimmer::processRecords(std::vector<neoReference>& recs, bool isPair, bool isReverse){
    for(neoReference& rec : recs){
        processOneRecord(rec);
    }
}

void AvgQualTrimmer::processRecords(std::vector<neoReference>& recs, int threadId, bool isPair, bool isReverse){
    for(neoReference& rec : recs){
        processOneRecord(rec);
    }
}
