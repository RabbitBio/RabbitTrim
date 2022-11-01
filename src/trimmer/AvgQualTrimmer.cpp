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
    int len = rec.lqual;
    char* seq = (char*)(rec.base +  rec.pseq);
    char* qual = (char*)(rec.base + rec.pqual);
    
    int total = 0;
    for(int i = 0; i < len; i++){
        int qual_val = seq[i] == 'N' ? 0 : qual[i] - phred;
        total += qual_val;
    }
    rec.lseq = total < qual * len ? 0 : len;
    rec.lqual = rec.lseq;

}


void AvgQualTrimmer::processRecords(std::vector<neoReference>& recs, bool isPair, bool isReverse){
    for(neoReference& rec : recs){
        processOneRecord(rec);
    }
}
