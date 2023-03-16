#include "trimmer/SeedClippingTrimmer.h"
#include <iostream>

// #if defined __SSE2__ && defined __AVX__ && defined __AVX2__ && defined TRIM_USE_VEC
#if defined __SSE2__  && defined TRIM_USE_VEC
#include <immintrin.h>
#endif
using namespace rabbit::trim;

SeedClippingTrimmer::SeedClippingTrimmer(double mismatch_, bool use_default_mismatch_, std::string seqA_, std::string seqB_, int minLen_, int minQual_, int window_, int phred_)
{
  mismatch = mismatch_;
  seqA = seqA_;
  seqB = seqB_;
  minLen = minLen_;
  minQual = minQual_;
  window = window_;
  phred = phred_;
  use_default_mismatch = use_default_mismatch_;

// #if defined __SSE2__ && defined __AVX__ && defined __AVX2__ && defined TRIM_USE_VEC
#if defined __SSE2__  && defined TRIM_USE_VEC
  std::cout << "using vectorisatino .. in Seed Clipping .. " << std::endl;
  seqA_str = new char[seqA.size() + 15];
  seqB_str = new char[seqB.size() + 15];
  for(int i = 0; i <  seqA.size(); i++)
  {
    seqA_str[i] = seqA[i];
  }
  for(int i = 0; i <  seqB.size(); i++)
  {
    seqB_str[i] = seqB[i];
  }

  // seqA(0, 2)
  index_1_1_str = new char[16];
  index_1_2_str = new char[16];
  index_1_3_str = new char[16];
  // seqB(0, 2)
  index_2_1_str = new char[16];
  index_2_2_str = new char[16];
  index_2_3_str = new char[16];
  // seqA(3, 5)
  index_3_1_str = new char[16];
  index_3_2_str = new char[16];
  index_3_3_str = new char[16];
  for(int i = 0; i < 16; i++)
  {
    index_1_1_str[i] = seqA[0];
    index_1_2_str[i] = seqA[1];
    index_1_3_str[i] = seqA[2];

    index_2_1_str[i] = seqB[0];
    index_2_2_str[i] = seqB[1];
    index_2_3_str[i] = seqB[2];

    index_3_1_str[i] = seqA[3];
    index_3_2_str[i] = seqA[4];
    index_3_3_str[i] = seqA[5];
  }
#endif
  
}

SeedClippingTrimmer::~SeedClippingTrimmer()
{
// #if defined __SSE2__ && defined __AVX__ && defined __AVX2__ && defined TRIM_USE_VEC
#if defined __SSE2__  && defined TRIM_USE_VEC
  delete [] seqA_str;
  delete [] seqB_str;

  delete [] index_1_1_str;
  delete [] index_1_2_str;
  delete [] index_1_3_str;

  delete [] index_2_1_str;
  delete [] index_2_2_str;
  delete [] index_2_3_str;

  delete [] index_3_1_str;
  delete [] index_3_2_str;
  delete [] index_3_3_str;
#endif
}

void SeedClippingTrimmer::processOneRecord(Reference& rec){}
void SeedClippingTrimmer::processRecords(std::vector<Reference>& recs, bool isPair, bool isReverse) {}

void SeedClippingTrimmer::processSingleRecord(neoReference& rec)
{
  rec.lorigin = 0; // 用lorigin区分是real_adapter(== 1) or tail_adapter(== 2)
  // low quality process 
  int rec_len = rec.lseq;
  char* rec_seq = (char*)(rec.base + rec.pseq);
  char* rec_qual = (char*)(rec.base + rec.pqual);
  int stop = minLen - 1;
  int i, j;
  for(i = rec_len - 1; i >= stop; )
  {
    // int qual = rec_seq[i] == 'N' ? 0 : (rec_qual[i] - phred);
    int qual = rec_qual[i] - phred;
    if(qual >= minQual)
    {
      for(j = i - 1; j > i - window; j--){
        if(j < 0 || (rec_qual[j] - phred < minQual)) break;
      }
      if(j == i - window) break;
      else i = j - 1;
    }
    else
      i--;
  }
  if(i < stop){
    rec.lseq = 0;
    rec.lqual = 0;
    return;
  } 
  rec.lseq = i + 1;
  rec.lqual = i + 1;
  rec_seq[rec.lseq] = '\0';
  
  // find seed
  std::set<int> seed;
  std::string tmp_index1 = seqA.substr(0, 3);
  std::string tmp_index2 = seqA.substr(3, 3);
  const char* adapter_index_1 = tmp_index1.c_str();
  const char* adapter_index_2 = tmp_index2.c_str();
#if defined __SSE2__ && defined TRIM_USE_VEC 
  rec_len = rec.lseq;
  int total = (rec_len - 2 + 13) / 14; // 14 is 16 - 2 if use xmm
  __m128i vindex1_1 = _mm_loadu_si128((__m128i_u*)(index_1_1_str)); 
  __m128i vindex1_2 = _mm_loadu_si128((__m128i_u*)(index_1_2_str)); 
  __m128i vindex1_3 = _mm_loadu_si128((__m128i_u*)(index_1_3_str)); 

  
  for(int i = 0; i < total; i++)
  {
      __m128i vrec = _mm_loadu_si128((__m128i_u*)(rec_seq + i * 14)); // SSE2 // need extra 16 Byte space 
      __m128i vres1_1  = _mm_cmpeq_epi8(vrec, vindex1_1);
      vrec = _mm_srli_si128(vrec, 0x1);
      __m128i vres1_2  = _mm_cmpeq_epi8(vrec, vindex1_2);
      vrec = _mm_srli_si128(vrec, 0x1);
      __m128i vres1_3  = _mm_cmpeq_epi8(vrec, vindex1_3);
      
      __m128i vres = _mm_and_si128(vres1_1, vres1_2);
      vres = _mm_and_si128(vres, vres1_3);

      int res = _mm_movemask_epi8(vres);
      if(__builtin_popcount(res))
      {
        for(int p = 0; p < 14; p++)
        {
          if((res >> p) & 1)
          {
            if(i * 14 + p <= rec_len - 3)
              seed.emplace(i * 14 + p);
          }
        }
        
      }
  }

  total = ((rec_len - 3) - 2 + 13 ) / 14; // 需要保证rec_len >= 3
  char* tmp_rec_seq = rec_seq + 3;
  __m128i vindex3_1 = _mm_loadu_si128((__m128i_u*)(index_3_1_str)); 
  __m128i vindex3_2 = _mm_loadu_si128((__m128i_u*)(index_3_2_str)); 
  __m128i vindex3_3 = _mm_loadu_si128((__m128i_u*)(index_3_3_str)); 
  for(int i = 0; i < total; i++)
  {
      __m128i vrec = _mm_loadu_si128((__m128i_u*)(tmp_rec_seq + i * 14)); // SSE2
      __m128i vres3_1  = _mm_cmpeq_epi8(vrec, vindex3_1);
      vrec = _mm_srli_si128(vrec, 0x1);
      __m128i vres3_2  = _mm_cmpeq_epi8(vrec, vindex3_2);
      vrec = _mm_srli_si128(vrec, 0x1);
      __m128i vres3_3  = _mm_cmpeq_epi8(vrec, vindex3_3);
      
      __m128i vres = _mm_and_si128(vres3_1, vres3_2);
      vres = _mm_and_si128(vres, vres3_3);

      int res = _mm_movemask_epi8(vres);
      if(__builtin_popcount(res))
      {
        for(int p = 0; p < 14; p++)
        {
          if((res >> p) & 1)
          {
            if(i * 14 + p <= rec_len - 6)
              seed.emplace(i * 14 + p);
          }
        }
        
      }
  }
  
#else
  char* indexloc = rec_seq;
  while (true)
  {
    indexloc = std::strstr(indexloc, adapter_index_1);
    if(indexloc != NULL){
      seed.emplace(indexloc - rec_seq);
      indexloc++;
    }
    else break;
  }

  // other adapter_index
  char* start = rec_seq + 3;
  indexloc = start;
  while (true)
  {
    indexloc = std::strstr(indexloc, adapter_index_2);
    if(indexloc != NULL){
      seed.emplace(indexloc - start);
      indexloc++;
    }
    else break;
  }
#endif

  //debug
  // std::cout << "seed is : \n";
  // for(auto iter = seed.begin(); iter != seed.end(); iter++){
  //   std::cout << *iter << " ";
  // }
  // std::cout << std::endl;

  // iterator seed set
  bool isFound;
  for(auto iter = seed.begin(); iter != seed.end(); iter++){
    int pos = *iter;
    int compLen = rec.lseq - pos; 
    compLen = compLen > seqA.size() ? seqA.size() : compLen;
    int maxMiss = use_default_mismatch ? ((compLen + 7) >> 3) : std::ceil(compLen * mismatch);
    char* start = rec_seq + pos;
    int tmpMiss = 0;
    isFound = true;

// #if defined __SSE2__ && defined __AVX__ && defined __AVX2__ && defined TRIM_USE_VEC
#if defined __SSE2__  && defined TRIM_USE_VEC
    
    int nums = (compLen + 15) / 16;
    int remain = compLen % 16;
    remain = remain == 0 ? 16 : remain;
    int mask = ((1 << remain) - 1);
    for(int i = 0; i < nums - 1; i++)
    {
      __m128i vrec = _mm_loadu_si128((__m128i_u*)(start + i * 16)); // SSE2
      __m128i vclip = _mm_loadu_si128((__m128i_u*)(seqA_str + i * 16));
      __m128i vres  = _mm_cmpeq_epi8(vrec, vclip);
      int res = _mm_movemask_epi8(vres); // SSE2
      int tmp_mis = 16 - __builtin_popcount(res);
      tmpMiss += tmp_mis;
    }
    __m128i vrec = _mm_loadu_si128((__m128i_u*)(start + (nums - 1) * 16));
    __m128i vclip = _mm_loadu_si128((__m128i_u*)(seqA_str + (nums - 1) * 16));
    __m128i vres  = _mm_cmpeq_epi8(vrec, vclip);
    int res = _mm_movemask_epi8(vres);
    res = res & mask;
    int tmp_mis = remain - __builtin_popcount(res);
    tmpMiss += tmp_mis;
    isFound = tmpMiss > maxMiss ? false : true;
    
#else
    for(int i = 0; i < compLen; i++){
      if(start[i] != seqA[i]){
        tmpMiss++;
        if(tmpMiss > maxMiss) 
        {
          isFound = false;
          break;
        }
      }
    }
#endif
    if(isFound){
      rec.lorigin = 1;
      if(pos >= minLen) rec.lseq = pos;
      else
      {
        rec.lseq = 0;
      if(pos <= DIMER_INSERT) rec.lorigin = 3; // 3 is real & dimer
      }
      rec.lqual = rec.lseq;
      return;
    }
  }

  if(rec_seq[rec.lseq - 2] == adapter_index_1[0] && rec_seq[rec.lseq - 1] == adapter_index_1[1])
  {
    rec.lorigin = 2;
    if(rec.lseq - 2 < minLen) rec.lseq = 0;
    else rec.lseq -= 2;
    rec.lqual = rec.lseq;
  }


}
void SeedClippingTrimmer::processPairRecord(neoReference& rec1, neoReference& rec2)
{
  rec1.lorigin = 0;
  rec2.lorigin = 0;

  int rec1_len = rec1.lseq;
  int rec2_len = rec2.lseq;
  char* rec1_seq = (char*)(rec1.base + rec1.pseq);
  char* rec2_seq = (char*)(rec2.base + rec2.pseq);
  char* rec1_qual = (char*)(rec1.base + rec1.pqual);
  char* rec2_qual = (char*)(rec2.base + rec2.pqual);

  // low quality trim
  int  rec_len = rec1_len < rec2_len ? rec1_len : rec2_len;
  int stop = minLen - 1;
  int i, j;
  int minQ = minQual + phred;
  for(i = rec_len - 1; i >= stop; )
  {
    if(rec1_qual[i] >= minQ && rec2_qual[i] >= minQ)
    {
      for(j = i - 1; j > i - window; j--)
      {
        if(j < 0 || rec1_qual[j] < minQ || rec2_qual[j] < minQ){
          break;
        }
      }
      if(j == i - window) break;
      else i = j - 1;
    }
    else i--;
  }

  if(i >= stop){
    rec_len = i + 1;
    rec1_seq[rec_len] = '\0';
    rec2_seq[rec_len] = '\0';
  }else{
    rec1.lseq = 0;
    rec1.lqual = 0;
    rec2.lseq = 0;
    rec2.lqual = 0;
    return;
  }

  // find seed
  std::set<int> seed;
  char* indexloc = rec1_seq;
  std::string tmp_index1 = seqA.substr(0, 3);
  std::string tmp_index2 = seqB.substr(0, 3);
  const char* adapter_index_1 = tmp_index1.c_str();
  const char* adapter_index_2 = tmp_index2.c_str();
  while(true)
  {
    indexloc = std::strstr(indexloc, adapter_index_1);
    if(indexloc != NULL) 
    {
      seed.emplace(indexloc - rec1_seq);
      indexloc++;
    }
    else break;
  }
  indexloc = rec2_seq;
  while(true)
  {
    indexloc = std::strstr(indexloc, adapter_index_2);
    if(indexloc != NULL) 
    {
      seed.emplace(indexloc - rec2_seq);
      indexloc++;
    }
    else break;
  }


  // iterator seed
  for(auto iter = seed.begin(); iter != seed.end(); iter++)
  {
    bool isFound = true;
    int pos = *iter;
    int compLen = rec_len - pos;
    compLen = compLen > seqA.size() ? seqA.size() : compLen;
    int maxMiss = (use_default_mismatch) ? ((compLen + 7) >> 3) : std::ceil(mismatch * compLen);
    char* rec1_start = rec1_seq + pos;
    char* rec2_start = rec2_seq + pos;
    int tmpMiss = 0;
#if defined __SSE2__  && defined TRIM_USE_VEC
    
    int nums = (compLen + 15) / 16;
    int remain = compLen % 16;
    remain = remain == 0 ? 16 : remain;
    int mask = ((1 << remain) - 1);
    for(int i = 0; i < nums - 1; i++)
    {
      __m128i vrec = _mm_loadu_si128((__m128i_u*)(rec1_start + i * 16)); // SSE2
      __m128i vclip = _mm_loadu_si128((__m128i_u*)(seqA_str + i * 16));
      __m128i vres  = _mm_cmpeq_epi8(vrec, vclip);
      int res = _mm_movemask_epi8(vres); // SSE2
      int tmp_mis = 16 - __builtin_popcount(res);
      tmpMiss += tmp_mis;
    }
    __m128i vrec = _mm_loadu_si128((__m128i_u*)(rec1_start + (nums - 1) * 16));
    __m128i vclip = _mm_loadu_si128((__m128i_u*)(seqA_str + (nums - 1) * 16));
    __m128i vres  = _mm_cmpeq_epi8(vrec, vclip);
    int res = _mm_movemask_epi8(vres);
    res = res & mask;
    int tmp_mis = remain - __builtin_popcount(res);
    tmpMiss += tmp_mis;
    isFound = tmpMiss > maxMiss ? false : true;
    
#else
    for(int i = 0; i < compLen; i++){
      if(rec1_start[i] != seqA[i]){
        tmpMiss++;
        if(tmpMiss > maxMiss){
          isFound = false;
          break;
        }
      }
    }

#endif

    if(isFound)
    { // rec1 has found
      tmpMiss = 0;
#if defined __SSE2__  && defined TRIM_USE_VEC
    
      int nums = (compLen + 15) / 16;
      int remain = compLen % 16;
      remain = remain == 0 ? 16 : remain;
      int mask = ((1 << remain) - 1);
      for(int i = 0; i < nums - 1; i++)
      {
        __m128i vrec = _mm_loadu_si128((__m128i_u*)(rec2_start + i * 16)); // SSE2
        __m128i vclip = _mm_loadu_si128((__m128i_u*)(seqB_str + i * 16));
        __m128i vres  = _mm_cmpeq_epi8(vrec, vclip);
        int res = _mm_movemask_epi8(vres); // SSE2
        int tmp_mis = 16 - __builtin_popcount(res);
        tmpMiss += tmp_mis;
      }
      __m128i vrec = _mm_loadu_si128((__m128i_u*)(rec2_start + (nums - 1) * 16));
      __m128i vclip = _mm_loadu_si128((__m128i_u*)(seqB_str + (nums - 1) * 16));
      __m128i vres  = _mm_cmpeq_epi8(vrec, vclip);
      int res = _mm_movemask_epi8(vres);
      res = res & mask;
      int tmp_mis = remain - __builtin_popcount(res);
      tmpMiss += tmp_mis;
      isFound = tmpMiss > maxMiss ? false : true;
    
#else
      for(int i = 0; i < compLen; i++)
      {
        if(rec2_start[i] != seqB[i])
        {
          tmpMiss++;
          if(tmpMiss > maxMiss)
          {
            isFound = false;
            break;
          }
        }
      }
#endif

      if(isFound)
      {
        // rec2 also found
        rec1.lorigin = 1;
        rec2.lorigin = 1;
        if(pos >= minLen){
          rec1.lseq = pos;
          rec1.lqual = pos;
          rec2.lseq = pos;
          rec2.lqual = pos;
        }else{
          if(pos <= DIMER_INSERT)
          {
            rec1.lorigin = 3;
            rec2.lorigin = 3;
          }
          rec1.lseq = 0;
          rec1.lqual = 0;
          rec2.lseq = 0;
          rec2.lqual = 0;
        }
        return;
      }
    }
  }

  // iter == seed.end() : all seed are not adapter
  ASSERT(rec_len >= 7);
  int tmpMiss = 0;
  if(rec1_seq[rec_len - 2] != seqA[0]) tmpMiss++;
  if(rec1_seq[rec_len - 1] != seqA[1]) tmpMiss++;
  if(rec2_seq[rec_len - 2] != seqB[0]) tmpMiss++;
  if(rec2_seq[rec_len - 1] != seqB[1]) tmpMiss++;
  if(!is_revcomp(rec1_seq[5], rec2_seq[rec_len - 8])) tmpMiss++;
  if(!is_revcomp(rec2_seq[5], rec1_seq[rec_len - 8])) tmpMiss++;
  if(tmpMiss <= 1)
  {
    rec1.lorigin = 2;
    rec2.lorigin = 2;
    rec_len -= 2;
    if(rec_len < minLen){
      rec_len = 0;
    }
  }
  else
  {
    tmpMiss = 0;
    if(rec1_seq[rec_len - 1] != seqA[0]) tmpMiss++;
    if(rec2_seq[rec_len - 1] != seqB[0]) tmpMiss++;
    if(!is_revcomp(rec1_seq[5], rec2_seq[rec_len - 7])) tmpMiss++;
    if(!is_revcomp(rec2_seq[5], rec1_seq[rec_len - 7])) tmpMiss++;
    if(!is_revcomp(rec1_seq[6], rec2_seq[rec_len - 8])) tmpMiss++;
    if(!is_revcomp(rec2_seq[6], rec1_seq[rec_len - 8])) tmpMiss++;
    if(tmpMiss <= 1){
      rec1.lorigin = 2;
      rec2.lorigin = 2;
      rec_len -= 1;
      if(rec_len < minLen) rec_len = 0;
    }
  }
  rec1.lseq  = rec_len;
  rec2.lseq  = rec_len;
  rec1.lqual = rec_len;
  rec2.lqual = rec_len;
}

void SeedClippingTrimmer::processOneRecord(neoReference& rec)
{

}
void SeedClippingTrimmer::processRecords(std::vector<neoReference>& recs, bool isPair, bool isReverse)
{
  if(isPair){
    int n = recs.size() / 2;
    ASSERT(recs.size() == n * 2);
    for(int i = 0; i < n; i++){
      neoReference& rec1 = recs[i];
      neoReference& rec2 = recs[i + n];
      processPairRecord(rec1, rec2);
    }
  }else{
    for(auto& rec : recs){
      processSingleRecord(rec);
    }
  }
}

void SeedClippingTrimmer::processRecords(std::vector<neoReference>& recs, int threadId, bool isPair, bool isReverse)
{
  if(isPair){
    int n = recs.size() / 2;
    ASSERT(recs.size() == n * 2);
    for(int i = 0; i < n; i++){
      neoReference& rec1 = recs[i];
      neoReference& rec2 = recs[i + n];
      processPairRecord(rec1, rec2);
    }
  }else{
    for(auto& rec : recs){
      processSingleRecord(rec);
    }
  }
}

