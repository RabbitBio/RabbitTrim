#include "trimmer/SeedClippingTrimmer.h"
#include <iostream>
using namespace rabbit::trim;

SeedClippingTrimmer::SeedClippingTrimmer(double mismatch_, std::string seqA_, std::string seqB_, int minLen_, int minQual_, int window_, int phred_)
{
  mismatch = mismatch_;
  seqA = seqA_;
  seqB = seqB_;
  minLen = minLen_;
  minQual = minQual_;
  window = window_;
  phred = phred_;
}

SeedClippingTrimmer::~SeedClippingTrimmer()
{

}

void SeedClippingTrimmer::processOneRecord(Reference& rec){}
void SeedClippingTrimmer::processRecords(std::vector<Reference>& recs, bool isPair, bool isReverse)
{
  for(auto& rec : recs){
    processOneRecord(rec);
  }
}



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

  // find seed
  std::set<int> seed;
  const char* adapter_index_1 = seqA.substr(0, 3).c_str();
  const char* adapter_index_2 = seqA.substr(3, 3).c_str();
  rec_seq[rec.lseq] = '\0';
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
    int maxMiss = std::ceil(compLen * mismatch);
    char* start = rec_seq + pos;
    int tmpMiss = 0;
    isFound = true;
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
    if(isFound){
      rec.lorigin = 1;
      if(pos >= minLen) rec.lseq = pos;
      else rec.lseq = 0;
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
  const char* adapter_index_1 = seqA.substr(0, 3).c_str();
  const char* adapter_index_2 = seqB.substr(0, 3).c_str();
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
    int maxMiss = std::ceil(mismatch * compLen);
    char* rec1_start = rec1_seq + pos;
    char* rec2_start = rec2_seq + pos;
    int tmpMiss = 0;
    for(int i = 0; i < compLen; i++){
      if(rec1_start[i] != seqA[i]){
        tmpMiss++;
        if(tmpMiss > maxMiss){
          isFound = false;
          break;
        }
      }
    }
    if(isFound)
    { // rec1 has found
      tmpMiss = 0;
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
      if(isFound)
      {
        rec1.lorigin = 1;
        rec2.lorigin = 1;
        if(pos >= minLen){
          rec1.lseq = pos;
          rec1.lqual = pos;
          rec2.lseq = pos;
          rec2.lqual = pos;
        }else{
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

