#include "trimmer/SlidingWindowTrimmer.h"
#include <vector>

using namespace rabbit::trim;

void SlidingWindowTrimmer::processOneRecord(Reference& rec){
    int len = rec.length;
    if(len < windowLength) {
        rec.length = 0;
        return;
    }
    int total = 0;
    std::vector<int> quals(len);
    for(int i = 0; i < windowLength; i++){
        int qual_val = (rec.seq.at(i + rec.headPos) == 'N') ? 0 : rec.quality.at(i + rec.headPos) - phred;
        quals[i] = qual_val;
        total += qual_val;
    }
    if(total < totalRequiredQuality) {
        rec.length = 0;
        return;
    }

    int lengthToKeep = len;
    for(int i = 0; i < len - windowLength; i++){
        int qual_val = (rec.seq.at(i+ rec.headPos + windowLength) == 'N') ? 0 : rec.quality.at(i + rec.headPos +windowLength) - phred;
        quals[i+windowLength] = qual_val;
        total = total - quals[i] + qual_val; 
        if(total < totalRequiredQuality){
            lengthToKeep = windowLength + i;
            break;
        }
    }

    // 滑动窗口结束后 从保留的最后一个位置开始再检查单个碱基的质量
    int k = lengthToKeep - 1;
    for(k; k >= 0; k--){
        if(quals[k] >= requiredQuality){
            lengthToKeep = k + 1;
            rec.length = lengthToKeep;
            return;
        }
    }

    rec.length = 0;
    
}
void SlidingWindowTrimmer::processRecords(std::vector<Reference>& recs, bool isPair, bool isReverse){
    for(Reference& rec : recs){
        processOneRecord(rec);
    }
}
void SlidingWindowTrimmer::processOneRecord(neoReference& rec){
    int len = rec.lseq;
    if(len < windowLength) {
        rec.lseq = 0;
        rec.lqual = 0;
        return;
    }
    int total = 0;
    std::vector<int> quals(len);
    char* rec_seq = (char*)(rec.base + rec.pseq);
    char* rec_qual = (char*)(rec.base + rec.pqual);
    for(int i = 0; i < windowLength; i++){
        int qual_val = (rec_seq[i] == 'N') ? 0 : rec_qual[i]- phred;
        quals[i] = qual_val;
        total += qual_val;
    }
    if(total < totalRequiredQuality) {
        rec.lseq = 0;
        rec.lqual = 0;
        return;
    }

    int lengthToKeep = len;
    for(int i = 0; i < len - windowLength; i++){
        int qual_val = (rec_seq[i + windowLength] == 'N') ? 0 : rec_qual[i + windowLength] - phred;
        quals[i+windowLength] = qual_val;
        total = total - quals[i] + qual_val; 
        if(total < totalRequiredQuality){
            lengthToKeep = windowLength + i;
            break;
        }
    }

    // 滑动窗口结束后 从保留的最后一个位置开始再检查单个碱基的质量
    int k = lengthToKeep - 1;
    for(k; k >= 0; k--){
        if(quals[k] >= requiredQuality){
            lengthToKeep = k + 1;
            rec.lseq = lengthToKeep;
            rec.lqual = lengthToKeep;
            return;
        }
    }

    rec.lseq = 0;
    rec.lqual = 0;
    
}
void SlidingWindowTrimmer::processRecords(std::vector<neoReference>& recs, bool isPair, bool isReverse){
    for(neoReference& rec : recs){
        processOneRecord(rec);
    }
}
void SlidingWindowTrimmer::processRecords(std::vector<neoReference>& recs, int threadId, bool isPair, bool isReverse){
    for(neoReference& rec : recs){
        processOneRecord(rec);
    }
}
