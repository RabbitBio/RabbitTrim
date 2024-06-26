#include "trimmer/TrailingTrimmer.h"

using namespace rabbit::trim;

void TrailingTrimmer::processOneRecord(Reference& rec){
    int len = rec.length;
    std::string seq = rec.seq;
    std::string quality = rec.quality;
    int cur_headPos = rec.headPos;
    for(int i = len - 1; i >= 0; i--){
        int qual_val = quality.at(i + cur_headPos) - phred;
        if(seq.at(i + cur_headPos) != 'N' && qual_val >= qual){
            // rec.seq = seq.substr(0, i+1);
            // rec.quality = quality.substr(0, i+1);
            rec.length = i + 1;
            return;
        }
    }
    rec.length = 0;
    // rec.seq = "";
    // rec.quality = "";
}


void TrailingTrimmer::processRecords(std::vector<Reference>& recs, bool isPair, bool isReverse){
    for(Reference& rec : recs){
        processOneRecord(rec);
    }
}

void TrailingTrimmer::processOneRecord(neoReference& rec){
    int len = rec.lseq;
    char* rec_qual = (char*)(rec.base + rec.pqual);
    char* rec_seq = (char*)(rec.base + rec.pseq);
    for(int i = len - 1; i >= 0; i--){
        int qual_val = rec_qual[i] - phred;
        if(rec_seq[i] != 'N' && qual_val >= qual){
            // rec.seq = seq.substr(0, i+1);
            // rec.quality = quality.substr(0, i+1);
            rec.lseq = i + 1;
            rec.lqual = i + 1;
            return;
        }
    }
    rec.lseq = 0;
    rec.lqual = 0;
    // rec.seq = "";
    // rec.quality = "";
}

void TrailingTrimmer::processRecords(std::vector<neoReference>& recs, bool isPair, bool isReverse){
    for(neoReference& rec : recs){
        processOneRecord(rec);
    }
}
void TrailingTrimmer::processRecords(std::vector<neoReference>& recs, int threadId, bool isPair, bool isReverse){
    for(neoReference& rec : recs){
        processOneRecord(rec);
    }
}
