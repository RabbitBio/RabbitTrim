#include "trimmer/IlluminaPrefixPair.h"

using namespace rabbit::trim;

IlluminaPrefixPair::IlluminaPrefixPair(std::string prefix1_, std::string prefix2_, rabbit::Logger& logger, int phred_, int minPrefix_, int seedMaxMiss_, int minPalindromeLikelihood_, bool palindromeKeepBoth_, int consumerNum_) : logger(logger){
  phred = phred_;
  minPrefix = minPrefix_;
  seedMaxMiss = seedMaxMiss_;
  seedMax = seedMaxMiss * 2;
  minPalindromeLikelihood = minPalindromeLikelihood_;
  palindromeKeepBoth = palindromeKeepBoth_;
  consumerNum = consumerNum_; 

  int prefix1Len = prefix1_.length();
  int prefix2Len = prefix2_.length();
  int minLength = prefix1Len < prefix2Len ? prefix1Len : prefix2Len;
  prefix1 = prefix1_.substr(0, minLength) ;
  prefix2 = prefix2_.substr(0, minLength) ;
  prefixLen = minLength;
  logger.infoln("Using PrefixPair: '" + prefix1 + "' and '" + prefix2 + "'"); // debug

  forwardPacks = new uint64[consumerNum_ * (rabbit::trim::MAX_READ_LENGTH + rabbit::trim::MAX_ADAPTER_LENGTH) ];
  reversePacks = new uint64[consumerNum_ * (rabbit::trim::MAX_READ_LENGTH + rabbit::trim::MAX_ADAPTER_LENGTH) ];
  // likelihoodArr = new float[consumerNum_ *  rabbit::trim::MAX_READ_LENGTH];
}

IlluminaPrefixPair::~IlluminaPrefixPair(){
  delete [] forwardPacks;
  delete [] reversePacks;
}

uint64 IlluminaPrefixPair::packCh(char ch, bool reverse){
  if(!reverse){
    switch (ch)
    {
      case 'A':
        return BASE_A;
        break;
      case 'C':
        return BASE_C;
        break;
      case 'G':
        return BASE_G;
        break;
      case 'T':
        return BASE_T;
        break;
    }

  }
  else{
    switch (ch)
    {
      case 'A':
        return BASE_T;
        break;
      case 'C':
        return BASE_G;
        break;
      case 'G':
        return BASE_C;
        break;
      case 'T':
        return BASE_A;
        break;
    }
  }
  return 0ULL;
}

//TODO 将Prefix 和 rec.seq 合并后再求packs, 和当前的方式比较看看哪一个效率更高
uint64* IlluminaPrefixPair::packSeqInternal(Reference& rec, bool reverse){
  int seqLen = rec.length;
  int len = prefixLen + seqLen;
  assert(len > 15);
  int cur_headPos = rec.headPos;
  uint64* out; // 用长度为16的窗口生产kmer
  uint64 pack = 0;
  out = new uint64[len - 15];
  if(!reverse){
    for(int i = 0; i < prefixLen; i++){
      uint64 tmp = packCh(prefix1.at(i), false);
      pack = (pack << 4) | tmp;
      if(i >= 15) out[i - 15] = pack; // 类似于长度为16的kmer
    }
    for(int i = 0; i < seqLen; i++) {
      uint64 tmp = packCh(rec.seq.at(i + cur_headPos), false);
      pack = (pack << 4) | tmp;
      if(i + prefixLen >= 15) out[i + prefixLen - 15] = pack; // 类似于长度为16的kmer
    }
  }
  else{
    for(int i = 0; i < prefixLen; i++){
      uint64 tmp = packCh(prefix2.at(i), true);
      pack = ((pack >> 4) | (tmp << 60));
      if(i >= 15) out[i - 15] = pack;
    }
    for(int i = 0; i < seqLen; i++) {
      uint64 tmp = packCh(rec.seq.at(i + cur_headPos), true);
      pack = (pack >> 4) | (tmp << 60); // 当rec是reverse read时，产生其反向互补链的pack
      if(i + prefixLen >= 15) out[i + prefixLen - 15] = pack; // 类似于长度为16的kmer
    }
  }
  return out;
}
uint64* IlluminaPrefixPair::packSeqInternal(neoReference& rec, bool reverse){
  int seqLen = rec.lseq;
  int len = prefixLen + seqLen;
  char* rec_seq = (char*)(rec.base + rec.pseq);
  assert(len > 15);
  uint64* out; // 用长度为16的窗口生产kmer
  uint64 pack = 0;
  out = new uint64[len - 15];
  if(!reverse){
    for(int i = 0; i < prefixLen; i++){
      uint64 tmp = packCh(prefix1.at(i), false);
      pack = (pack << 4) | tmp;
      if(i >= 15) out[i - 15] = pack; // 类似于长度为16的kmer
    }
    for(int i = 0; i < seqLen; i++) {
      uint64 tmp = packCh(rec_seq[i], false);
      pack = (pack << 4) | tmp;
      if(i + prefixLen >= 15) out[i + prefixLen - 15] = pack; // 类似于长度为16的kmer
    }
  }
  else{
    for(int i = 0; i < prefixLen; i++){
      uint64 tmp = packCh(prefix2.at(i), true);
      pack = ((pack >> 4) | (tmp << 60));
      if(i >= 15) out[i - 15] = pack;
    }
    for(int i = 0; i < seqLen; i++) {
      uint64 tmp = packCh(rec_seq[i], true);
      pack = (pack >> 4) | (tmp << 60); // 当rec是reverse read时，产生其反向互补链的pack
      if(i + prefixLen >= 15) out[i + prefixLen - 15] = pack; // 类似于长度为16的kmer
    }
  }
  return out;
}

uint64* IlluminaPrefixPair::packSeqInternalForward(neoReference& rec, int threadId){
  int seqLen = rec.lseq;
  int len = prefixLen + seqLen;
  char* rec_seq = (char*)(rec.base + rec.pseq);
  assert(len > 15);
  uint64 pack = 0;
  uint64* out = forwardPacks + (rabbit::trim::MAX_READ_LENGTH + rabbit::trim::MAX_ADAPTER_LENGTH) * threadId;
  for(int i = 0; i < prefixLen; i++){
    uint64 tmp = (1 << ((prefix1[i] >> 1) & 7)) & 15;
    pack = (pack << 4) | tmp;
    out[i] = pack; 
  }
  for(int i = 0; i < seqLen; i++) {
    uint64 tmp = (1 << ((rec_seq[i] >> 1) & 7)) & 15;
    pack = (pack << 4) | tmp;
    out[i + prefixLen] = pack;
  }
  return out + 15;
}

uint64* IlluminaPrefixPair::packSeqInternalReverse(neoReference& rec, int threadId){
  int seqLen = rec.lseq;
  int len = prefixLen + seqLen;
  char* rec_seq = (char*)(rec.base + rec.pseq);
  assert(len > 15);
  uint64 pack = 0;
  uint64* out = reversePacks + (rabbit::trim::MAX_READ_LENGTH + rabbit::trim::MAX_ADAPTER_LENGTH) * threadId;
  for(int i = 0; i < prefixLen; i++){
    uint64 forwardCoding = ((1 << ((prefix2[i] >> 1) & 7)) & 15); 
    uint64 tmp = ((forwardCoding >> 2) | (forwardCoding << 2)) & 15;
    pack = ((pack >> 4) | (tmp << 60));
    out[i] = pack; 
  }
  for(int i = 0; i < seqLen; i++) {
    uint64 forwardCoding = ((1 << ((rec_seq[i] >> 1) & 7)) & 15); 
    uint64 tmp = ((forwardCoding >> 2) | (forwardCoding << 2)) & 15;
    pack = ((pack >> 4) | (tmp << 60));
    out[i + prefixLen] = pack;
  }
  return out + 15;
}


float IlluminaPrefixPair::calculatePalindromeDifferenceQuality(Reference& rec1, Reference& rec2, int overlap, int skip1, int skip2){
  // compute quality
  int len1 = rec1.length;
  int cur_headPos1 = rec1.headPos;
  int len2 = rec2.length;
  int cur_headPos2 = rec2.headPos;

  float likelihood;
  float totalLikelihood = 0;

  for(int i = 0; i < overlap; i++){
    int offset1 = i + skip1;
    int offset2 = skip2 + overlap - i - 1;

    // TODO 向量化
    char ch1 = offset1 < prefixLen ? prefix1.at(offset1) : rec1.seq.at(cur_headPos1 + offset1 - prefixLen);
    char ch2 = offset2 < prefixLen ? prefix2.at(offset2) : rec2.seq.at(cur_headPos2 + offset2 - prefixLen);


    int qual1 = offset1 < prefixLen ? 100 : (rec1.quality.at(cur_headPos1 + offset1 - prefixLen) - phred);
    int qual2 = offset2 < prefixLen ? 100 : (rec2.quality.at(cur_headPos2 + offset2 - prefixLen) - phred);
    int minQual = qual1 < qual2 ? qual1 : qual2;

    float s = ((ch1 >> 1) & 3) == (((ch2 >> 1) & 3) ^ 2) ? LOG10_4 :  -minQual / 10.0f; // XOR 2表示取碱基的互补碱基

    likelihood = (ch1 == 'N' || ch2 == 'N') ? 0 : s;
    totalLikelihood += likelihood;
  }
  return totalLikelihood;
}

float IlluminaPrefixPair::calculatePalindromeDifferenceQuality(neoReference& rec1, neoReference& rec2, int overlap, int skip1, int skip2){
  // compute quality
  int len1 = rec1.lseq;
  int len2 = rec2.lseq;
  char* rec1_seq = (char*)(rec1.base + rec1.pseq);
  char* rec2_seq = (char*)(rec2.base + rec2.pseq);
  char* rec1_qual = (char*)(rec1.base + rec1.pqual);
  char* rec2_qual = (char*)(rec2.base + rec2.pqual);

  float likelihood;
  float totalLikelihood = 0;

  for(int i = 0; i < overlap; i++){
    int offset1 = i + skip1;
    int offset2 = skip2 + overlap - i - 1;

    // TODO 向量化
    char ch1 = offset1 < prefixLen ? prefix1.at(offset1) : rec1_seq[offset1 - prefixLen];
    char ch2 = offset2 < prefixLen ? prefix2.at(offset2) : rec2_seq[offset2 - prefixLen];


    int qual1 = offset1 < prefixLen ? 100 : (rec1_qual[offset1 - prefixLen] - phred);
    int qual2 = offset2 < prefixLen ? 100 : (rec2_qual[offset2 - prefixLen] - phred);
    int minQual = qual1 < qual2 ? qual1 : qual2;

    float s = ((ch1 >> 1) & 3) == (((ch2 >> 1) & 3) ^ 2) ? LOG10_4 :  -minQual / 10.0f; // XOR 2表示取碱基的互补碱基

    likelihood = (ch1 == 'N' || ch2 == 'N') ? 0 : s;
    totalLikelihood += likelihood;
  }
  return totalLikelihood;
}

float IlluminaPrefixPair::calculatePalindromeDifferenceQuality(neoReference& rec1, neoReference& rec2, int overlap, int skip1, int skip2, int threadId){
  // compute quality
  int len1 = rec1.lseq;
  int len2 = rec2.lseq;
  char* rec1_seq = (char*)(rec1.base + rec1.pseq);
  char* rec2_seq = (char*)(rec2.base + rec2.pseq);
  char* rec1_qual = (char*)(rec1.base + rec1.pqual);
  char* rec2_qual = (char*)(rec2.base + rec2.pqual);

  float likelihood;
  float totalLikelihood = 0;
  // float* likelihood_arr = likelihoodArr + threadId * rabbit::trim::MAX_READ_LENGTH;
  // int cnt = 0;
  
  int offset1 = skip1;
  int offset2 = skip2 + overlap - 1;
  int minCmpLen, maxCmpLen;
  // if(prefixLen - skip1 < skip2 + overlap - prefixLen){
    minCmpLen = std::max(prefixLen - skip1, 0);
    maxCmpLen = std::min(skip2 + overlap - prefixLen, overlap);
  // }
  // else{
  //   maxCmpLen = prefixLen - skip1;
  //   minCmpLen = skip2 + overlap - prefixLen;
  // }

  for(int i = 0; i < minCmpLen; i++)
  {
    char ch1 = prefix1[offset1];
    char ch2 = rec2_seq[offset2 - prefixLen];
    int minQual = rec2_qual[offset2 - prefixLen] - phred;

    float s = ((ch1 >> 1) & 3) == (((ch2 >> 1) & 3) ^ 2) ? LOG10_4 :  -minQual / 10.0f; // XOR 2表示取碱基的互补碱基

    likelihood = (ch1 == 'N' || ch2 == 'N') ? 0.0f : s;
    totalLikelihood += likelihood;
    // likelihood_arr[cnt++] = likelihood;
    offset1++;
    offset2--;
  }
  
  if(prefixLen - skip1 == minCmpLen)
  {
    for(int i = 0; i < maxCmpLen - minCmpLen; i++)
    {
      char ch1 = rec1_seq[offset1 - prefixLen];
      char ch2 = rec2_seq[offset2 - prefixLen];
      int qual1 = rec1_qual[offset1 - prefixLen] - phred;
      int qual2 = rec1_qual[offset2 - prefixLen] - phred;
      int minQual = qual1 < qual2 ? qual1 : qual2;

      float s = ((ch1 >> 1) & 3) == (((ch2 >> 1) & 3) ^ 2) ? LOG10_4 :  -minQual / 10.0f; // XOR 2表示取碱基的互补碱基
      likelihood = (ch1 == 'N' || ch2 == 'N') ? 0 : s;
      totalLikelihood += likelihood;
      // likelihood_arr[cnt++] = likelihood;
      offset1++;
      offset2--;
    }
    
    for(int i = 0; i < overlap - maxCmpLen; i++)
    {
      char ch1 = rec1_seq[offset1 - prefixLen];
      char ch2 = prefix2[offset2];
      int minQual = rec1_qual[offset1 - prefixLen] - phred;
      float s = ((ch1 >> 1) & 3) == (((ch2 >> 1) & 3) ^ 2) ? LOG10_4 :  -minQual / 10.0f; // XOR 2表示取碱基的互补碱基

      likelihood = (ch1 == 'N' || ch2 == 'N') ? 0 : s;
      totalLikelihood += likelihood;
      // likelihood_arr[cnt++] = likelihood;
      offset1++;
      offset2--;
    }

  }
  else
  {
    logger.errorln("not impossiblbe");
    // is is possible? No, wmk think
    for(int i = 0; i < maxCmpLen - minCmpLen; i++)
    {
      char ch1 = prefix1[offset1]; 
      char ch2 = prefix2[offset2];
      int minQual = 100;

      float s = ((ch1 >> 1) & 3) == (((ch2 >> 1) & 3) ^ 2) ? LOG10_4 :  -minQual / 10.0f; // XOR 2表示取碱基的互补碱基
      likelihood = (ch1 == 'N' || ch2 == 'N') ? 0 : s;
      totalLikelihood += likelihood;
      offset1++;
      offset2--;
    }
    
    for(int i = 0; i < overlap - maxCmpLen; i++)
    {
      char ch1 = rec1_seq[offset1 - prefixLen];
      char ch2 = prefix2[offset2];
      int minQual = rec1_qual[offset1 - prefixLen] - phred;
      float s = ((ch1 >> 1) & 3) == (((ch2 >> 1) & 3) ^ 2) ? LOG10_4 :  -minQual / 10.0f; // XOR 2表示取碱基的互补碱基

      likelihood = (ch1 == 'N' || ch2 == 'N') ? 0 : s;
      totalLikelihood += likelihood;
      offset1++;
      offset2--;
    }
  }
  // for(int i = 0; i < cnt; i++)
  // {
  //   totalLikelihood += likelihood_arr[i];
  // }
  return totalLikelihood;
}


int IlluminaPrefixPair::palindromeReadsCompare(Reference& rec1, Reference& rec2){
  // prefix + rec.seq 的长度为16的kmer数组 ， 注意：对于reverse read得到的是其反向互补的kmer 
  uint64* pack1 = packSeqInternal(rec1, false);
  uint64* pack2 = packSeqInternal(rec2, true);

  int testIndex = 0;
  int refIndex = prefixLen; // 参考序列从不包含prefix的第一个kmer开始

  int rec1Len = rec1.length;
  int rec2Len = rec2.length;

  int pack1Len = rec1Len + prefixLen - 15;
  int pack2Len = rec2Len + prefixLen - 15;


  // if(pack1.length <= refIndex)
  // if(rec1Len <= 15 || rec2Len <= 15) return std::INT_MAX; 
  if(rec1Len <= 15 || rec2Len <= 15) return 1 << 30; // why?

  int count = 0;
  int seedSkip = prefixLen - 16; 
  if(seedSkip > 0) {
    testIndex = seedSkip;
    count = seedSkip;
  }

  uint64 ref1, ref2;

  int seqlen1 = rec1Len + prefixLen;
  int seqlen2 = rec2Len + prefixLen; 

  int maxCount = (seqlen1 > seqlen2 ? seqlen1 : seqlen2) - 15 - minPrefix;

  // count = testIndex + refIndex - prefixLen

  for(int i = count; i < maxCount; i++){
    ref1 = pack1[refIndex];
    ref2 = pack2[refIndex]; 
    if((testIndex < pack2Len && __builtin_popcountll(ref1 ^ pack2[testIndex]) <= seedMax) || (testIndex < pack1Len && __builtin_popcountll(ref2 ^ pack1[testIndex]) <= seedMax)){
      int totalOverLap = count + prefixLen + 16;
      int skip1 = 0;
      int skip2 = 0;

      if(totalOverLap > seqlen2) skip1 = totalOverLap - seqlen2;
      if(totalOverLap > seqlen1) skip2 = totalOverLap - seqlen1;

      int actualOverlap = totalOverLap - skip1 - skip2;

      float palindromeLikelihood = calculatePalindromeDifferenceQuality(rec1, rec2, actualOverlap, skip1, skip2);

      if(palindromeLikelihood >= minPalindromeLikelihood){
        assert(totalOverLap - 2 * prefixLen);
        return totalOverLap - 2 * prefixLen ; // 返回的是不包括prefix的两条序列overlap的长度 即序列的真实长度 NOTE:返回值可能大于序
      }
    }
    count++;
    if((count & 1) == 0 && refIndex + 1 < pack1Len  && refIndex + 1 < pack2Len) refIndex++; // count是偶数移动refIndex 奇数移动testIndex
    else testIndex++;
  }
  return 1 << 30;
}

int IlluminaPrefixPair::palindromeReadsCompare(neoReference& rec1, neoReference& rec2){
  // prefix + rec.seq 的长度为16的kmer数组 ， 注意：对于reverse read得到的是其反向互补的kmer 
  uint64* pack1 = packSeqInternal(rec1, false);
  uint64* pack2 = packSeqInternal(rec2, true);

  int testIndex = 0;
  int refIndex = prefixLen; // 参考序列从不包含prefix的第一个kmer开始

  int rec1Len = rec1.lseq;
  int rec2Len = rec2.lseq;

  int pack1Len = rec1Len + prefixLen - 15;
  int pack2Len = rec2Len + prefixLen - 15;


  // if(pack1.length <= refIndex)
  // if(rec1Len <= 15 || rec2Len <= 15) return std::INT_MAX; 
  if(rec1Len <= 15 || rec2Len <= 15) return 1 << 30; // why?

  int count = 0;
  int seedSkip = prefixLen - 16; 
  if(seedSkip > 0) {
    testIndex = seedSkip;
    count = seedSkip;
  }

  uint64 ref1, ref2;

  int seqlen1 = rec1Len + prefixLen;
  int seqlen2 = rec2Len + prefixLen; 

  int maxCount = (seqlen1 > seqlen2 ? seqlen1 : seqlen2) - 15 - minPrefix;

  // count = testIndex + refIndex - prefixLen

  for(int i = count; i < maxCount; i++){
    ref1 = pack1[refIndex];
    ref2 = pack2[refIndex]; 
    if((testIndex < pack2Len && __builtin_popcountll(ref1 ^ pack2[testIndex]) <= seedMax) || (testIndex < pack1Len && __builtin_popcountll(ref2 ^ pack1[testIndex]) <= seedMax)){
      int totalOverLap = count + prefixLen + 16;
      int skip1 = 0;
      int skip2 = 0;

      if(totalOverLap > seqlen2) skip1 = totalOverLap - seqlen2;
      if(totalOverLap > seqlen1) skip2 = totalOverLap - seqlen1;

      int actualOverlap = totalOverLap - skip1 - skip2;

      float palindromeLikelihood = calculatePalindromeDifferenceQuality(rec1, rec2, actualOverlap, skip1, skip2);

      if(palindromeLikelihood >= minPalindromeLikelihood){
        assert(totalOverLap - 2 * prefixLen);
        return totalOverLap - 2 * prefixLen ; // 返回的是不包括prefix的两条序列overlap的长度 即序列的真实长度 NOTE:返回值可能大于序
      }
    }
    count++;
    if((count & 1) == 0 && refIndex + 1 < pack1Len  && refIndex + 1 < pack2Len) refIndex++; // count是偶数移动refIndex 奇数移动testIndex
    else testIndex++;
  }
  return 1 << 30;
}

int IlluminaPrefixPair::palindromeReadsCompare(neoReference& rec1, neoReference& rec2, int threadId){
  int rec1Len = rec1.lseq;
  int rec2Len = rec2.lseq;
  if(rec1Len <= 15 || rec2Len <= 15) return 1 << 30; // why?

  // prefix + rec.seq 的长度为16的kmer数组 ， 注意：对于reverse read得到的是其反向互补的kmer 
  uint64* pack1 = packSeqInternalForward(rec1, threadId); 
  uint64* pack2 = packSeqInternalReverse(rec2, threadId); 

  int testIndex = 0;
  int refIndex = prefixLen; // 参考序列从不包含prefix的第一个kmer开始

  int pack1Len = rec1Len + prefixLen - 15;
  int pack2Len = rec2Len + prefixLen - 15;

  int count = 0;
  int seedSkip = prefixLen - 16; 
  if(seedSkip > 0) {
    testIndex = seedSkip;
    count = seedSkip;
  }

  uint64 ref1, ref2;

  int seqlen1 = rec1Len + prefixLen;
  int seqlen2 = rec2Len + prefixLen; 

  int maxCount = (seqlen1 > seqlen2 ? seqlen1 : seqlen2) - 15 - minPrefix;

  // count = testIndex + refIndex - prefixLen

  for(int i = count; i < maxCount; i++){
    ref1 = pack1[refIndex];
    ref2 = pack2[refIndex]; 
    if((testIndex < pack2Len && __builtin_popcountll(ref1 ^ pack2[testIndex]) <= seedMax) || (testIndex < pack1Len && __builtin_popcountll(ref2 ^ pack1[testIndex]) <= seedMax)){
      int totalOverLap = count + prefixLen + 16;
      int skip1 = 0;
      int skip2 = 0;

      if(totalOverLap > seqlen2) skip1 = totalOverLap - seqlen2;
      if(totalOverLap > seqlen1) skip2 = totalOverLap - seqlen1;

      int actualOverlap = totalOverLap - skip1 - skip2;

      float palindromeLikelihood = calculatePalindromeDifferenceQuality(rec1, rec2, actualOverlap, skip1, skip2, threadId);

      if(palindromeLikelihood >= minPalindromeLikelihood){
        assert(totalOverLap - 2 * prefixLen);
        return totalOverLap - 2 * prefixLen ; // 返回的是不包括prefix的两条序列overlap的长度 即序列的真实长度 NOTE:返回值可能大于序
      }
    }
    count++;
    if((count & 1) == 0 && refIndex + 1 < pack1Len  && refIndex + 1 < pack2Len) refIndex++; // count是偶数移动refIndex 奇数移动testIndex
    else testIndex++;
  }
  return 1 << 30;
}
