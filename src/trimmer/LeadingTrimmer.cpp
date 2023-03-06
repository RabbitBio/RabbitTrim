#include "trimmer/LeadingTrimmer.h"
#include <string>

using namespace rabbit::trim;
    
void LeadingTrimmer::processOneRecord(Reference& rec){
    std::string quality = rec.quality;
    std::string seq = rec.seq;
    int len = rec.length;
    int cur_headPos = rec.headPos;
    for(int i = 0; i < len; i++){
        // N 的质量分数当做0
        char ch = seq.at(i + cur_headPos);
        int qual_val = quality.at(i + cur_headPos) - phred;
        if(ch != 'N' && qual_val >= qual){
            // rec.seq = seq.substr(i);
            // rec.quality = quality.substr(i);
            rec.headPos += i;
            rec.length = len - i;
            return;
        }
        
    }
    // The entire sequence quality score is not met
    // rec.quality = "";
    // rec.seq = "";
    rec.headPos += len;
    rec.length = 0;
}

void LeadingTrimmer::processRecords(std::vector<Reference>& recs, bool isPair, bool isReverse){
    for(Reference& rec : recs){
        processOneRecord(rec);
    }
}

void LeadingTrimmer::processOneRecord(neoReference& rec){
    char* rec_qual = (char*)(rec.base + rec.pqual);
    char* rec_seq = (char*)(rec.base + rec.pseq);
    int len = rec.lseq;
    for(int i = 0; i < len; i++){
        // N 的质量分数当做0
        char ch = rec_seq[i];
        int qual_val = rec_qual[i] - phred;
        if(ch != 'N' && qual_val >= qual){
            // rec.seq = seq.substr(i);
            // rec.quality = quality.substr(i);
            rec.pseq += i;
            rec.pqual += i;
            rec.lseq = len - i;
            rec.lqual = rec.lseq;
            return;
        }
        
    }
    // The entire sequence quality score is not met
    // rec.quality = "";
    // rec.seq = "";
    rec.pseq += len;
    rec.pqual += len;
    rec.lseq = 0;
    rec.lqual = 0;
}

void LeadingTrimmer::processRecords(std::vector<neoReference>& recs, bool isPair, bool isReverse){
    for(neoReference& rec : recs){
        processOneRecord(rec);
    }
}
void LeadingTrimmer::processRecords(std::vector<neoReference>& recs, int threadId, bool isPair, bool isReverse){
    for(neoReference& rec : recs){
        processOneRecord(rec);
    }
}
