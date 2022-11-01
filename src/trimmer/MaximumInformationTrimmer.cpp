#include "trimmer/MaximumInformationTrimmer.h"

using namespace rabbit::trim;

double MaximumInformationTrimmer::calcNormalization(double* arr, int arrLength, int margin){
    double maxVal = arr[0];
    for(int i = 1; i < arrLength; i++){
        double val = std::abs(arr[i]);
        maxVal = val > maxVal ? val : maxVal;
    }
    
    return std::numeric_limits<long long>::max() / (maxVal * margin);
}

int64* MaximumInformationTrimmer::normalize(double* arr, int arrLength, double ratio){
    int64* out = new int64[arrLength] ;
    for(int i = 0; i < arrLength; i++){
        out[i] = (int64)(arr[i] * ratio);
    }
    return out;
}

MaximumInformationTrimmer::MaximumInformationTrimmer(int parLength_, float strictness_){
    parLength = parLength_;
    strictness = strictness_;
    
    // Calculate the score for each length in advance
    lengthScoreTmp = new double[LONGEST_READ];
    for(int i = 0; i < LONGEST_READ; i++){
        // length threshold factor
        double pow1 = std::exp(parLength - i - 1);
        double unique = std::log(1.0 / (1.0 + pow1));
        
        // coverage factor
        double coverage = std::log(i+1) * (1 - strictness);

        lengthScoreTmp[i] = unique + coverage;
        
    }
    
    qualProbTmp = new double[MAXQUAL + 1];
    for(int i = 0; i < MAXQUAL + 1; i++){
        // [0,60] -> [0.109, 1) before log
        // [0,60] -> [-2.2164, 0) * strictness (after log)
        // error rate factor
        qualProbTmp[i] = std::log(1 - std::pow(0.1, (0.5+i)/10.0)) * strictness;
    }
    
    double normRatio = std::max(calcNormalization(lengthScoreTmp, LONGEST_READ, LONGEST_READ * 2), calcNormalization(qualProbTmp, MAXQUAL+1, LONGEST_READ * 2));

    lengthScore = normalize(lengthScoreTmp, LONGEST_READ, normRatio);
    qualProb = normalize(qualProbTmp, 1 + MAXQUAL, normRatio);
    
}
void MaximumInformationTrimmer::processOneRecord(Reference& rec){
    int len = rec.length;
    // compute quality 
    int cur_headPos = rec.headPos;
    int* quals = new int[len];
    for(int i = 0; i < len; i++) {
        int qual_val = rec.seq.at(i + cur_headPos) == 'N' ? 0 : rec.quality.at(i + cur_headPos) - phred;
        quals[i] = qual_val;
    }

    // Accumulated quality score
    int64 accumQuality = 0;
    int64 maxScore = std::numeric_limits<long long>::min();

    int maxScorePosition = 0;
    for(int i = 0; i < len; i++){
        int q = quals[i];
        q = q < 0 ? 0 : q;
        q = q > MAXQUAL ? MAXQUAL : q;
        accumQuality += qualProb[q];
        
        int64 score = lengthScore[i] + accumQuality;
        // maxScore = score > maxScore ? score : maxScore;
        if(score >= maxScore){
            maxScore = score;
            maxScorePosition = i; 
        }
        
    }
    // maxScore == 0 ? // TODO 
    rec.length = maxScorePosition + 1;
    
}

void MaximumInformationTrimmer::processRecords(std::vector<Reference>& recs, bool isPair, bool isReverse){
    for(Reference& rec : recs){
        processOneRecord(rec);
    }
}

void MaximumInformationTrimmer::processOneRecord(neoReference& rec){
    int len = rec.lseq;
    char* rec_seq = (char*)(rec.base + rec.pseq);
    char* rec_qual = (char*)(rec.base + rec.pqual);
    // compute quality 
    int* quals = new int[len];
    for(int i = 0; i < len; i++) {
        int qual_val = rec_seq[i] == 'N' ? 0 : rec_qual[i] - phred;
        quals[i] = qual_val;
    }

    // Accumulated quality score
    int64 accumQuality = 0;
    int64 maxScore = std::numeric_limits<long long>::min();

    int maxScorePosition = 0;
    for(int i = 0; i < len; i++){
        int q = quals[i];
        q = q < 0 ? 0 : q;
        q = q > MAXQUAL ? MAXQUAL : q;
        accumQuality += qualProb[q];
        
        int64 score = lengthScore[i] + accumQuality;
        // maxScore = score > maxScore ? score : maxScore;
        if(score >= maxScore){
            maxScore = score;
            maxScorePosition = i; 
        }
        
    }
    // maxScore == 0 ? // TODO 
    rec.lseq = maxScorePosition + 1;
    rec.lqual = rec.lseq;
    
}

void MaximumInformationTrimmer::processRecords(std::vector<neoReference>& recs, bool isPair, bool isReverse){
    for(neoReference& rec : recs){
        processOneRecord(rec);
    }
}
