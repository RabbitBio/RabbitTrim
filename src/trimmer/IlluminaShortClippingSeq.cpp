#include "trimmer/IlluminaShortClippingSeq.h"

#define likely(x) __builtin_expect(!!(x),1)
#define unlikely(x) __builtin_expect(!!(x),0)

using namespace rabbit::trim;

inline uint64 IlluminaShortClippingSeq::calcSingleMask(int length){
    uint64 mask = 0xffffffffffffffffULL; 
    if(length < 16)
        mask <<= (16 - length) * 4;
    return mask;
}

float IlluminaShortClippingSeq::calculateMaximumRange(float* vals, int valsLen){
    // 找到vals子区间元素和的最大值
    float sum = 0;
    float max = vals[0];
    for(int i = 0; i < valsLen; i++ ){
        if(sum < 0) sum = vals[i];
        else sum += vals[i];
        if(sum > max) max = sum;
    }
    return max;
}
// TODO
// inline float calculateMaximumRange(float* vals, int valsLen){
//     float total = 0;
//     std::vector<float> merges;
//     for(int i = 0; i < valsLen; i++){
//         int val = vals[i];
//         if((total > 0 && val < 0) || (total < 0 && val > 0)){
//             merges.emplace_back(total);
//             total = val;
//         }
//         else{
//             total += val;
//         }
//     }
//     merges.emplace_back(total);
//     
//     bool scanAgain = true;
//     while(!merges.empty() && scanAgain){
//         scanAgain = false;
//         auto mergeIter = merges.begin();
//         for(mergeIter; mergeIter != merges.end(); mergeIter++){
//             float val = *mergeIter;  
//             
//         }
//     }
// }


IlluminaShortClippingSeq::IlluminaShortClippingSeq(rabbit::Logger& logger_, int phred_, std::string seq_, int seedMaxMiss_, int minSequenceLikelihood_, int minSequenceOverlap_, int consumerNum_ ): logger(logger_){
    phred = phred_; 
    seedMaxMiss = seedMaxMiss_;
    seedMax = seedMaxMiss * 2;
    minSequenceLikelihood = minSequenceLikelihood_;
    minSequenceOverlap = minSequenceOverlap_;
    consumerNum = consumerNum_;

    seq = seq_;
    logger.infoln("Using Short Clipping Sequence: '" + seq + "'");
    seqLen = seq.length();
    // calcSingleMask 
    assert(seqLen < 16);
    mask = 0xffffffffffffffffULL;
    mask <<= ((16 - seqLen) * 4);
    
    // packSeqExternal clipPack
    pack = new uint64[seqLen];
    uint64 pack_ = 0ULL;
    for(int i = 0; i < 15 + seqLen ; i++){
        int tmp = 0;
        if(i < seqLen)
            tmp = packCh(seq.at(i), false);
        pack_ = (pack_ << 4) | tmp; 
        if(i >= 15) pack[i - 15] = pack_;
    }
    
    // for each thread, create one rec pack
    recPacks = new uint64[(rabbit::trim::MAX_READ_LENGTH) * consumerNum_];
}

IlluminaShortClippingSeq::~IlluminaShortClippingSeq(){
    delete [] pack;
    delete [] recPacks;
}


int IlluminaShortClippingSeq::readsSeqCompare(Reference& rec){
    uint64* packRec = packSeqExternal(rec);
    int packRecLen = rec.length; 
    uint64* packClip = pack;
    int packClipLen = seqLen;

    int packRecMax = packRecLen - minSequenceOverlap;
    int packClipMax = packClipLen - minSequenceOverlap; // 是不是应该+1呀？ 不然可能为0
    
    
    std::set<int> offsetSet;
    // offset 默认升序排列 

    for(int i = 0; i < packRecMax; i++){
        uint64 comboMask = calcSingleMask(packRecLen - i) & mask;
        for(int j = 0; j < packClipMax; j++){
            // TODO 为什么不重新计算Clip的mask
            int diff = __builtin_popcountll((packRec[i] ^ packClip[j]) & comboMask);
            if(diff <= seedMax){
                int offset = i - j;
                offsetSet.emplace(offset);
            }
        }
    }
    
    for(std::set<int>::iterator iter = offsetSet.begin(); iter != offsetSet.end(); iter++){
        // Iterate through offsetSet from smallest to largest
        int offset = *iter;
        int recCompLength = offset > 0 ? rec.length - offset : rec.length;
        int clipCompLength = offset < 0 ? seqLen + offset : seqLen;
        int compLength = recCompLength < clipCompLength ? recCompLength : clipCompLength;

        // debug  compLength 是不是一定大于minSequenceOverlap
        assert(compLength > minSequenceOverlap);
        if(compLength > minSequenceOverlap){
            float seqLikelihood = calculateDifferenceQuality(rec, compLength, offset);
            if(seqLikelihood >= minSequenceLikelihood) return offset;
        }
        
    }
    // return std::INT_MAX;
    return 1 << 30;
}

int IlluminaShortClippingSeq::readsSeqCompare(neoReference& rec, int threadId){
    uint64* packRec = packSeqExternal(rec, threadId);
    int packRecLen = rec.lseq; 
    uint64* packClip = pack;
    int packClipLen = seqLen;

    int packRecMax = packRecLen - minSequenceOverlap;
    int packClipMax = packClipLen - minSequenceOverlap; // 是不是应该+1呀？ 不然可能为0
    
    
    // offset 从可能的最小值到最大值依次选择 保持i j 的相对关系寻找 直到第一个offset满足 (效果不好...
    // robin_hood::unordered_set<int> offsetSet;
    // int rv = 0 - (packClipMax - 1);
    // int maxRelativeVal = packRecMax - 1;
    // 
    // while(rv <= maxRelativeVal)
    // {
    //   int i = 0;
    //   int j = i - rv;
    //   while(i < packRecMax && j < packClipMax)
    //   {
    //     uint64 comboMask = calcSingleMask(packRecLen - i) & mask;
    //     int diff = __builtin_popcountll((packRec[i] ^ packClip[j]) & comboMask);
    //     if(diff <= seedMax)
    //     {
    //       offsetSet.emplace(rv);
    //       break;
    //     }
    //     else
    //     {
    //       i++;
    //       j++;
    //     }
    //     
    //   }
    //   rv++;
    // }


    int minOffset = 0 - (packClipMax - 1);
    int maxOffset = packRecMax - 1;
    int offsetNums = maxOffset - minOffset + 1;
    int cntLen = (offsetNums + 63) / 64;
    uint64* cnt = new uint64[cntLen];
    for(int i = 0; i < cntLen; i++)
    {
      cnt[i] = 0ULL;
    }

    for(int i = 0; i < packRecMax; i++){
        uint64 comboMask = calcSingleMask(packRecLen - i) & mask;
        for(int j = 0; j < packClipMax/4*4; j+=4){
            int diff1 = __builtin_popcountll((packRec[i] ^ packClip[j]) & comboMask);
            int diff2 = __builtin_popcountll((packRec[i] ^ packClip[j+1]) & comboMask);
            int diff3 = __builtin_popcountll((packRec[i] ^ packClip[j+2]) & comboMask);
            int diff4 = __builtin_popcountll((packRec[i] ^ packClip[j+3]) & comboMask);
            if(unlikely(diff1 <= seedMax)){
              int offset = i - j;
              // offsetSet.emplace(offset);
              cnt[((offset - minOffset) >> 6)] |= 1ULL << (63 - ((offset - minOffset) % 64));
            }
            if(unlikely(diff2 <= seedMax)){
              int offset = i - j - 1;
              // offsetSet.emplace(offset);
              cnt[((offset - minOffset) >> 6)] |= 1ULL << (63 - ((offset - minOffset) % 64));
            }
            if(unlikely(diff3 <= seedMax)){
              int offset = i - j - 2;
              // offsetSet.emplace(offset);
              cnt[((offset - minOffset) >> 6)] |= 1ULL << (63 - ((offset - minOffset) % 64));
            }
            if(unlikely(diff4 <= seedMax)){
              int offset = i - j - 3;
              // offsetSet.emplace(offset);
              cnt[((offset - minOffset) >> 6)] |= 1ULL << (63 - ((offset - minOffset) % 64));
            }
        }
        for(int j = packClipMax / 4 * 4; j < packClipMax; j++)
        {
            int diff = __builtin_popcountll((packRec[i] ^ packClip[j]) & comboMask);
            if(unlikely(diff <= seedMax)){
              int offset = i - j;
              // offsetSet.emplace(offset);
              cnt[((offset - minOffset) >> 6)] |= 1ULL << (63 - ((offset - minOffset) % 64));
            }

        }
    }

    // std::set<int> offsetSet;
    // // offset 默认升序排列 
    // for(int i = 0; i < packRecMax; i++){
    //     uint64 comboMask = calcSingleMask(packRecLen - i) & mask;
    //     for(int j = 0; j < packClipMax; j++){
    //         // TODO 为什么不重新计算Clip的mask
    //         int diff = __builtin_popcountll((packRec[i] ^ packClip[j]) & comboMask);
    //         if(diff <= seedMax){
    //             int offset = i - j;
    //             total++;
    //             // if(offsets.count(offset)){
    //             //   cnt++;
    //             //   continue;
    //             // }
    //             // offsets.emplace(offset);
    //             offsetSet.emplace(offset);
    //         }
    //     }
    // }
    
    // for(auto iter = offsetSet.begin(); iter != offsetSet.end(); iter++){
    for(int i = 0; i < cntLen; i++){
        // Iterate through offsetSet from smallest to largest
        for(int p = 0; p < 64; p++){
          if((cnt[i] >> (63 - p))&1){
            int offset = minOffset + i * 64 + p;
            int recCompLength = offset > 0 ? rec.lseq - offset : rec.lseq;
            int clipCompLength = offset < 0 ? seqLen + offset : seqLen;
            int compLength = recCompLength < clipCompLength ? recCompLength : clipCompLength;

            // debug  compLength 是不是一定大于minSequenceOverlap
            assert(compLength > minSequenceOverlap);
            if(compLength > minSequenceOverlap){
              float seqLikelihood = calculateDifferenceQuality(rec, compLength, offset);
              if(seqLikelihood >= minSequenceLikelihood) return offset;
            }
          }
        }
    }
    // return std::INT_MAX;
    return 1 << 30;
}

int IlluminaShortClippingSeq::readsSeqCompare(neoReference& rec){
    uint64* packRec = packSeqExternal(rec);
    int packRecLen = rec.lseq; 
    uint64* packClip = pack;
    int packClipLen = seqLen;

    int packRecMax = packRecLen - minSequenceOverlap;
    int packClipMax = packClipLen - minSequenceOverlap; // 是不是应该+1呀？ 不然可能为0
    
    
    // offset 从可能的最小值到最大值依次选择 保持i j 的相对关系寻找 直到第一个offset满足 效果不好...
    // robin_hood::unordered_set<int> offsetSet;
    // int rv = 0 - (packClipMax - 1);
    // int maxRelativeVal = packRecMax - 1;
    // 
    // while(rv <= maxRelativeVal)
    // {
    //   int i = 0;
    //   int j = i - rv;
    //   while(i < packRecMax && j < packClipMax)
    //   {
    //     uint64 comboMask = calcSingleMask(packRecLen - i) & mask;
    //     int diff = __builtin_popcountll((packRec[i] ^ packClip[j]) & comboMask);
    //     if(diff <= seedMax)
    //     {
    //       offsetSet.emplace(rv);
    //       break;
    //     }
    //     else
    //     {
    //       i++;
    //       j++;
    //     }
    //     
    //   }
    //   rv++;
    // }


    int minOffset = 0 - (packClipMax - 1);
    int maxOffset = packRecMax - 1;
    int offsetNums = maxOffset - minOffset + 1;
    int cntLen = (offsetNums + 63) / 64;
    uint64* cnt = new uint64[cntLen];
    for(int i = 0; i < cntLen; i++)
    {
      cnt[i] = 0ULL;
    }

    for(int i = 0; i < packRecMax; i++){
        uint64 comboMask = calcSingleMask(packRecLen - i) & mask;
        for(int j = 0; j < packClipMax; j++){
            // TODO 为什么不重新计算Clip的mask
            int diff = __builtin_popcountll((packRec[i] ^ packClip[j]) & comboMask);
            if(diff <= seedMax){
              int offset = i - j;
              // offsetSet.emplace(offset);
              cnt[((offset - minOffset) >> 6)] |= 1ULL << (63 - ((offset - minOffset) % 64));
            }
        }
    }

    // std::set<int> offsetSet;
    // // offset 默认升序排列 
    // for(int i = 0; i < packRecMax; i++){
    //     uint64 comboMask = calcSingleMask(packRecLen - i) & mask;
    //     for(int j = 0; j < packClipMax; j++){
    //         // TODO 为什么不重新计算Clip的mask
    //         int diff = __builtin_popcountll((packRec[i] ^ packClip[j]) & comboMask);
    //         if(diff <= seedMax){
    //             int offset = i - j;
    //             total++;
    //             // if(offsets.count(offset)){
    //             //   cnt++;
    //             //   continue;
    //             // }
    //             // offsets.emplace(offset);
    //             offsetSet.emplace(offset);
    //         }
    //     }
    // }
    
    // for(auto iter = offsetSet.begin(); iter != offsetSet.end(); iter++){
    for(int i = 0; i < cntLen; i++){
        // Iterate through offsetSet from smallest to largest
        for(int p = 0; p < 64; p++){
          if((cnt[i] >> (63 - p))&1){
            int offset = minOffset + i * 64 + p;
            int recCompLength = offset > 0 ? rec.lseq - offset : rec.lseq;
            int clipCompLength = offset < 0 ? seqLen + offset : seqLen;
            int compLength = recCompLength < clipCompLength ? recCompLength : clipCompLength;

            // debug  compLength 是不是一定大于minSequenceOverlap
            assert(compLength > minSequenceOverlap);
            if(compLength > minSequenceOverlap){
              float seqLikelihood = calculateDifferenceQuality(rec, compLength, offset);
              if(seqLikelihood >= minSequenceLikelihood) return offset;
            }
          }
        }
    }
    // return std::INT_MAX;
    return 1 << 30;
}

// 用长度为16的窗口，从序列第一个位置开始，到最后一个位置结束，当剩余序列长度<16时，补充0
uint64* IlluminaShortClippingSeq::packSeqExternal(Reference& rec){ 
  int len = rec.length;
  int cur_headPos = rec.headPos;
  uint64* out = new uint64[len];

  uint64 pack = 0ULL;

  for(int i = 0; i < len + 15; i++){
    uint64 tmp = 0;
    if(i < len)
      tmp = packCh(rec.seq.at(cur_headPos + i), false);
    pack = (pack << 4) | tmp;
    if(i >= 15) out[i - 15] = pack;
  }
  return out;
}
uint64* IlluminaShortClippingSeq::packSeqExternal(neoReference& rec){ 
  int len = rec.lseq;
  char* rec_seq = (char*)(rec.base +  rec.pseq);
  uint64* out = new uint64[len];
  uint64 pack = 0ULL;

  for(int i = 0; i < len + 15; i++){
    uint64 tmp = 0;
    if(i < len)
      tmp = packCh(rec_seq[i], false);
    // tmp = 1 << ((rec_seq[i] >> 1) & 3); // no effect..
    pack = (pack << 4) | tmp;
    if(i >= 15) out[i - 15] = pack;
  }
  return out;
}

// using memory pool
uint64* IlluminaShortClippingSeq::packSeqExternal(neoReference& rec, int threadId){ 
  int len = rec.lseq;
  char* rec_seq = (char*)(rec.base +  rec.pseq);
  uint64* out = recPacks + threadId * MAX_READ_LENGTH;
  uint64 pack = 0ULL;

  for(int i = 0; i < len + 15; i++){
    uint64 tmp = 0;
    if(i < len)
      tmp = packCh(rec_seq[i], false);
    pack = (pack << 4) | tmp;
    if(i >= 15) out[i - 15] = pack;
  }
  return out;
}


int IlluminaShortClippingSeq::packCh(char ch, bool reverse){
  if(!reverse){
    switch (ch)
    {
      case 'A':
        return BASE_A;
      case 'C':
        return BASE_C;
      case 'G':
        return BASE_G;
      case 'T':
        return BASE_T;
    }

  }
  else{
    switch (ch)
    {
      case 'A':
        return BASE_T;
      case 'C':
        return BASE_G;
      case 'G':
        return BASE_C;
      case 'T':
        return BASE_A;
    }
  }
  return 0;
}

float IlluminaShortClippingSeq::calculateDifferenceQuality(Reference& rec, int overlap, int recOffset){
  // getQualityAsInteger
  int len = rec.length;
  int cur_headPos = rec.headPos;

  // Location to start comparison
  int recPos = recOffset > 0 ? recOffset : 0;
  int clipPos = recOffset < 0 ? -recOffset : 0;

  float* likelihood = new float[overlap];
  for(int i = 0; i < overlap; i++){
    char ch1 = rec.seq.at(cur_headPos + recPos);
    char ch2 = seq.at(clipPos);

    int qual_val = rec.quality.at(cur_headPos + recPos) - phred;
        float s = ((ch1 >> 1) & 3) == ((ch2 >> 1) & 3) ? LOG10_4 : -qual_val / 10.0f;
        likelihood[i] = ((ch1 == 'N' || ch2 == 'N') ? 0 : s);
        recPos++;
        clipPos++;
    }
    // calculateMaximumRange()
    float l = calculateMaximumRange(likelihood, overlap);
    return l;
}

float IlluminaShortClippingSeq::calculateDifferenceQuality(neoReference& rec, int overlap, int recOffset){
    // getQualityAsInteger
    int len = rec.lseq;
    char* rec_seq = (char*)(rec.base + rec.pseq);
    char* rec_qual = (char*)(rec.base + rec.pqual);

    // Location to start comparison
    int recPos = recOffset > 0 ? recOffset : 0;
    int clipPos = recOffset < 0 ? -recOffset : 0;

    float* likelihood = new float[overlap];
    for(int i = 0; i < overlap; i++){
        char ch1 = rec_seq[recPos];
        char ch2 = seq.at(clipPos);
        
        int qual_val = rec_qual[recPos] - phred;
        float s = ((ch1 >> 1) & 3) == ((ch2 >> 1) & 3) ? LOG10_4 : -qual_val / 10.0f;
        likelihood[i] = ((ch1 == 'N' || ch2 == 'N') ? 0 : s);
        recPos++;
        clipPos++;
    }
    // calculateMaximumRange()
    float l = calculateMaximumRange(likelihood, overlap);
    return l;
}
