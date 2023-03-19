#include "trimmer/MaximumInformationTrimmer.h"

#if defined __SSE2__ && defined __SSE4_1__ && defined TRIM_USE_VEC
#include <immintrin.h>
#endif

using namespace rabbit::trim;


MaximumInformationTrimmer::MaximumInformationTrimmer(int parLength_, float strictness_, int phred_, int consumerNum_){
    parLength = parLength_;
    strictness = strictness_;
    phred = phred_;

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

    qualsTotal = new char[consumerNum_ * rabbit::trim::MAX_READ_LENGTH];

    // std::cout << "lengthScore : " << std::endl;
    // for(int i = 0; i < LONGEST_READ; i++)
    // {
    //   std::cout << lengthScore[i] << std::endl;
    // }
    // std::cout << "qualProb : " << std::endl;
    // for(int i = 0; i < 1 + MAXQUAL; i++)
    // {
    //   std::cout << qualProb[i] << std::endl;
    // }

#if defined __SSE2__ && defined __SSE4_1__ && defined TRIM_USE_VEC
    all_N = new char[16];
    phred_arr = new char[16];
    max_qual = new char[16];
    for(int i = 0; i < 16; i++)
    {
      all_N[i] = 'N';
      phred_arr[i] = phred_;
      max_qual[i] = 60;
    }

#endif
    
}
MaximumInformationTrimmer::~MaximumInformationTrimmer()
{
  delete [] lengthScoreTmp;
  delete [] qualProbTmp;
  delete [] lengthScore;
  delete [] qualProb;
  delete [] qualsTotal;
#if defined __SSE2__ && defined __SSE4_1__ && defined TRIM_USE_VEC
  delete [] all_N;
  delete [] phred_arr;
  delete [] max_qual;
#endif
  
}

void MaximumInformationTrimmer::processOneRecord(neoReference& rec){}

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


void MaximumInformationTrimmer::processOneRecord(neoReference& rec, int threadId){
    int len = rec.lseq;
    char* rec_seq = (char*)(rec.base + rec.pseq);
    char* rec_qual = (char*)(rec.base + rec.pqual);
    // compute quality 
    char* quals = qualsTotal + threadId * rabbit::trim::MAX_READ_LENGTH;

#if defined __SSE2__ && defined __SSE4_1__ && defined TRIM_USE_VEC
    int nums = (len + 15) / 16;
    for(int i = 0; i < nums; i++)
    {
      
      __m128i vphred = _mm_loadu_si128((__m128i_u*)(phred_arr));
      __m128i vmaxq = _mm_loadu_si128((__m128i_u*)max_qual);
      __m128i vallN = _mm_loadu_si128((__m128i_u*)(all_N));
      __m128i vrec = _mm_loadu_si128((__m128i_u*)(rec_seq + i * 16));
      __m128i vqual = _mm_loadu_si128((__m128i_u*)(rec_qual + i * 16));

      __m128i vres1 = _mm_sub_epi8(vqual, vphred);
      __m128i vzero = _mm_setzero_si128();
      vres1 = _mm_max_epi8(vres1, vzero);
      vres1 = _mm_min_epi8(vres1, vmaxq);

      __m128i vres2 = _mm_cmpeq_epi8(vrec, vallN);

      __m128i vres = _mm_blendv_epi8(vres1, vzero, vres2); // SSE 4.1
      _mm_storeu_si128((__m128i_u*)(quals + i * 16), vres);
    }
    
#else
    for(int i = 0; i < len; i++) {
        quals[i] = rec_seq[i] == 'N' ? 0 : rec_qual[i] - phred;
    }
#endif

    // Accumulated quality score
    int64 accumQuality = 0;
    int64 maxScore = std::numeric_limits<long long>::min();

    int maxScorePosition = 0;
    for(int i = 0; i < len; i++){
        int q = quals[i];
#if defined __SSE2__ && defined __SSE4_1__ && defined TRIM_USE_VEC
#else
        q = q < 0 ? 0 : q;
        q = q > MAXQUAL ? MAXQUAL : q;
#endif
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
    if(maxScorePosition < 1 || maxScore == 0) rec.lseq = 0;
    rec.lqual = rec.lseq;
    
}

void MaximumInformationTrimmer::processRecords(std::vector<neoReference>& recs, bool isPair, bool isReverse){}

void MaximumInformationTrimmer::processRecords(std::vector<Reference>& recs, bool isPair, bool isReverse){
    for(Reference& rec : recs){
        processOneRecord(rec);
    }
}

void MaximumInformationTrimmer::processRecords(std::vector<neoReference>& recs, int threadId, bool isPair, bool isReverse){
    for(neoReference& rec : recs){
        processOneRecord(rec, threadId);
    }
}

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
