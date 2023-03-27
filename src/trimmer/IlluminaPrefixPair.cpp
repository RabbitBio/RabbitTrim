#include "trimmer/IlluminaPrefixPair.h"
#if defined __SSE2__ && defined __SSSE3__ && defined __AVX__ && defined __AVX2__ && defined TRIM_USE_VEC
#include <immintrin.h>
#endif

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

#if defined __SSE2__ && defined __SSSE3__ && defined __AVX__ && defined __AVX2__ && defined TRIM_USE_VEC
  all_N = new char[16];
  all_4 = new char[16];
  all_6 = new char[16];
  phred_arr = new char[16];
  index_arr = new char[16];
  for(int i = 0; i < 16; i++)
  {
    all_N[i] = 'N';
    all_4[i] =  4;
    all_6[i] =  6;
    phred_arr[i] = phred;
    index_arr[i] = 15 - i;
  }
  divide_arr = new float[8];
  awards = new float[8];
  for(int i = 0; i < 8; i++)
  {
    divide_arr[i] = -10.0f;
    awards[i] = LOG10_4;
  }
  likelihoodTotal = new float[ rabbit::trim::MAX_READ_LENGTH * consumerNum_];
  
#endif
}

IlluminaPrefixPair::~IlluminaPrefixPair(){
  delete [] forwardPacks;
  delete [] reversePacks;
#if defined __SSE2__ && defined __SSSE3__ && defined __AVX__ && defined __AVX2__ && defined TRIM_USE_VEC
  delete [] all_N;
  delete [] phred_arr;
  delete [] index_arr;
  delete [] divide_arr;
  delete [] awards;
  delete [] likelihoodTotal;
  delete [] all_4;
  delete [] all_6;
  
#endif
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

    float s = ((ch1 >> 1) & 3) == (((ch2 >> 1) & 3) ^ 2) ? LOG10_4 :  -minQual / 10; // XOR 2表示取碱基的互补碱基

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

    float s = ((ch1 >> 1) & 3) == (((ch2 >> 1) & 3) ^ 2) ? LOG10_4 :  -minQual / 10; // XOR 2表示取碱基的互补碱基

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
  // int cnt = 0;
  
  int offset1 = skip1;
  int offset2 = skip2 + overlap - 1;
  int minCmpLen, maxCmpLen;
  // if(prefixLen - skip1 < skip2 + overlap - prefixLen){ //  must go here
    // minCmpLen = std::max(prefixLen - skip1, 0);
    // maxCmpLen = std::min(skip2 + overlap - prefixLen, overlap);
  minCmpLen = prefixLen - skip1 > 0 ? prefixLen - skip1 : 0;
  maxCmpLen = skip2 - prefixLen < 0 ? skip2 + overlap - prefixLen : overlap;
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

    float s = ((ch1 >> 1) & 3) == (((ch2 >> 1) & 3) ^ 2) ? LOG10_4 :  -minQual / 10; // XOR 2表示取碱基的互补碱基

    likelihood = (ch1 == 'N' || ch2 == 'N') ? 0.0f : s;
    totalLikelihood += likelihood;
    offset1++;
    offset2--;
  }
  
#if defined __SSE2__ && defined __SSSE3__ && defined __AVX__ && defined __AVX2__ && defined TRIM_USE_VEC
  // 16B extra space required for front and rear each
  int tmp_overlap = maxCmpLen - minCmpLen;
  int tmp_pos_1 = skip1 + minCmpLen - prefixLen;
  int tmp_pos_2 = skip2 + overlap - 1 - minCmpLen - prefixLen - 15;
  char* tmp_rec_seq_1 = rec1_seq + tmp_pos_1;
  char* tmp_rec_seq_2 = rec2_seq + tmp_pos_2;
  char* tmp_rec_qual_1 = rec1_qual + tmp_pos_1;
  char* tmp_rec_qual_2 = rec2_qual + tmp_pos_2;
  float* tmp_likelihood = likelihoodTotal + threadId * rabbit::trim::MAX_READ_LENGTH;
  int nums = (tmp_overlap + 15) / 16;
  
  for(int i = 0; i < nums; i++)
  {
    __m128i vindex  = _mm_loadu_si128((__m128i_u*)(index_arr));
    __m128i vqual1  = _mm_loadu_si128((__m128i_u*)(tmp_rec_qual_1 + i * 16));
    __m128i vqual2  = _mm_loadu_si128((__m128i_u*)(tmp_rec_qual_2 - i * 16));

    vqual2 = _mm_shuffle_epi8(vqual2, vindex); // ssse3
    __m128i vqual = _mm_min_epu8(vqual1, vqual2); // sse2
    __m128i vphred  = _mm_loadu_si128((__m128i_u*)phred_arr);
    __m128i vscore  = _mm_sub_epi8(vqual, vphred);
    __m128i vscore_sr8 = _mm_srli_si128(vscore, 0x8);
    __m256i vscore_1 = _mm256_cvtepi8_epi32(vscore);
    __m256i vscore_2 = _mm256_cvtepi8_epi32(vscore_sr8);
    __m256 vscore_1_ps = _mm256_cvtepi32_ps(vscore_1);
    __m256 vscore_2_ps = _mm256_cvtepi32_ps(vscore_2);
    __m256 vdivide = _mm256_loadu_ps(divide_arr);
    vscore_1_ps = _mm256_div_ps(vscore_1_ps, vdivide);
    vscore_2_ps = _mm256_div_ps(vscore_2_ps, vdivide);

    vscore_1_ps = _mm256_ceil_ps(vscore_1_ps);
    vscore_2_ps = _mm256_ceil_ps(vscore_2_ps);

    __m128i vrec1_origin  = _mm_loadu_si128((__m128i_u*)(tmp_rec_seq_1 + i * 16));
    __m128i vrec2_origin  = _mm_loadu_si128((__m128i_u*)(tmp_rec_seq_2 - i * 16));

    __m128i vrec2_reverse = _mm_shuffle_epi8(vrec2_origin, vindex); // ssse3
    __m128i vall4 = _mm_loadu_si128((__m128i_u*)(all_4));
    __m128i vall6 = _mm_loadu_si128((__m128i_u*)(all_6));
    
    __m128i vrec1 = _mm_and_si128(vrec1_origin,  vall6);
    __m128i vrec2 = _mm_and_si128(vrec2_reverse, vall6);
    vrec2 = _mm_xor_si128(vrec2, vall4);

    __m128i vres  = _mm_cmpeq_epi8(vrec1, vrec2);
    __m128i vres_sr8 = _mm_srli_si128(vres, 0x8);
    __m256i vlow = _mm256_cvtepi8_epi32(vres);
    __m256i vhigh = _mm256_cvtepi8_epi32(vres_sr8);
    __m256 vlow_ps = _mm256_cvtepi32_ps(vlow);
    __m256 vhigh_ps = _mm256_cvtepi32_ps(vhigh);
    __m256 vawards = _mm256_loadu_ps(awards);
    vscore_1_ps = _mm256_blendv_ps(vscore_1_ps, vawards, vlow_ps);
    vscore_2_ps = _mm256_blendv_ps(vscore_2_ps, vawards, vhigh_ps);


    // judge N
    __m128i vallN = _mm_loadu_si128((__m128i_u*)(all_N));
    __m128i vres1  = _mm_cmpeq_epi8(vrec1_origin, vallN);
    __m128i vres1_sr8 = _mm_srli_si128(vres1, 0x8);
    __m256i vlow1 = _mm256_cvtepi8_epi32(vres1);
    __m256i vhigh1 = _mm256_cvtepi8_epi32(vres1_sr8);
    __m256 vlow1_ps = _mm256_cvtepi32_ps(vlow1);
    __m256 vhigh1_ps = _mm256_cvtepi32_ps(vhigh1);
    __m256 vscore_0 = _mm256_setzero_ps();
    vscore_1_ps = _mm256_blendv_ps(vscore_1_ps, vscore_0, vlow1_ps);
    vscore_2_ps = _mm256_blendv_ps(vscore_2_ps, vscore_0, vhigh1_ps);

    __m128i vres2  = _mm_cmpeq_epi8(vrec2_reverse, vallN);
    __m128i vres2_sr8 = _mm_srli_si128(vres2, 0x8);
    __m256i vlow2 = _mm256_cvtepi8_epi32(vres2);
    __m256i vhigh2 = _mm256_cvtepi8_epi32(vres2_sr8);
    __m256 vlow2_ps = _mm256_cvtepi32_ps(vlow2);
    __m256 vhigh2_ps = _mm256_cvtepi32_ps(vhigh2);
    vscore_1_ps = _mm256_blendv_ps(vscore_1_ps, vscore_0, vlow2_ps);
    vscore_2_ps = _mm256_blendv_ps(vscore_2_ps, vscore_0, vhigh2_ps);

    _mm256_storeu_ps(tmp_likelihood + i * 16,      vscore_1_ps);
    _mm256_storeu_ps(tmp_likelihood + i * 16 + 8,  vscore_2_ps);
  }

  for(int i = 0; i < maxCmpLen - minCmpLen; i++)
    totalLikelihood += tmp_likelihood[i]; 
  
  
  offset1 += tmp_overlap;
  offset2 -= tmp_overlap;

#else
  
  // if(prefixLen - skip1 == minCmpLen)
  // { // must go here
    for(int i = 0; i < maxCmpLen - minCmpLen; i++)
    {
      char ch1 = rec1_seq[offset1 - prefixLen];
      char ch2 = rec2_seq[offset2 - prefixLen];
      int qual1 = rec1_qual[offset1 - prefixLen] - phred;
      int qual2 = rec2_qual[offset2 - prefixLen] - phred;
      int minQual = qual1 < qual2 ? qual1 : qual2;

      // float s = ((ch1 >> 1) & 3) == (((ch2 >> 1) & 3) ^ 2) ? LOG10_4 :  -minQual / 10; // XOR 2表示取碱基的互补碱基
      float s = ((ch1) & 6) == (((ch2) & 6) ^ 4) ? LOG10_4 :  -minQual / 10; // XOR 2表示取碱基的互补碱基
      likelihood = (ch1 == 'N' || ch2 == 'N') ? 0 : s;
      totalLikelihood += likelihood;
      offset1++;
      offset2--;
    }
#endif

    for(int i = 0; i < overlap - maxCmpLen; i++)
    {
      char ch1 = rec1_seq[offset1 - prefixLen];
      char ch2 = prefix2[offset2];
      int minQual = rec1_qual[offset1 - prefixLen] - phred;
      float s = ((ch1 >> 1) & 3) == (((ch2 >> 1) & 3) ^ 2) ? LOG10_4 :  -minQual / 10; // XOR 2表示取碱基的互补碱基

      likelihood = (ch1 == 'N' || ch2 == 'N') ? 0 : s;
      totalLikelihood += likelihood;
      offset1++;
      offset2--;
    }

  // }
  // else
  // {
  //   logger.errorln("impossiblbe");
  //   // is is possible? No, wmk think
  //   for(int i = 0; i < maxCmpLen - minCmpLen; i++)
  //   {
  //     char ch1 = prefix1[offset1]; 
  //     char ch2 = prefix2[offset2];
  //     int minQual = 100;

  //     float s = ((ch1 >> 1) & 3) == (((ch2 >> 1) & 3) ^ 2) ? LOG10_4 :  -minQual / 10; // XOR 2表示取碱基的互补碱基
  //     likelihood = (ch1 == 'N' || ch2 == 'N') ? 0 : s;
  //     totalLikelihood += likelihood;
  //     offset1++;
  //     offset2--;
  //   }
  //   
  //   for(int i = 0; i < overlap - maxCmpLen; i++)
  //   {
  //     char ch1 = rec1_seq[offset1 - prefixLen];
  //     char ch2 = prefix2[offset2];
  //     int minQual = rec1_qual[offset1 - prefixLen] - phred;
  //     float s = ((ch1 >> 1) & 3) == (((ch2 >> 1) & 3) ^ 2) ? LOG10_4 :  -minQual / 10; // XOR 2表示取碱基的互补碱基

  //     likelihood = (ch1 == 'N' || ch2 == 'N') ? 0 : s;
  //     totalLikelihood += likelihood;
  //     offset1++;
  //     offset2--;
  //   }
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
