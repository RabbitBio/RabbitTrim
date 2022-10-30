#include "IlluminaShortClippingSeq.h"

using namespace rabbit::trim;

uint64 IlluminaShortClippingSeq::calcSingleMask(int length){
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


IlluminaShortClippingSeq::IlluminaShortClippingSeq(rabbit::Logger& logger_, int phred_, std::string seq_, int seedMaxMiss_, int minSequenceLikelihood_, int minSequenceOverlap_ ): logger(logger_){
    phred = phred_; 
    seedMaxMiss = seedMaxMiss_;
    seedMax = seedMaxMiss * 2;
    minSequenceLikelihood = minSequenceLikelihood_;
    minSequenceOverlap = minSequenceOverlap_;
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
}

IlluminaShortClippingSeq::~IlluminaShortClippingSeq(){
    delete [] pack;
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

// 用长度为16的窗口，从序列第一个位置开始，到最后一个位置结束，当剩余序列长度<16时，补充0
uint64* IlluminaShortClippingSeq::packSeqExternal(Reference& rec){ 
    int len = rec.length;
    int cur_headPos = rec.headPos;
    uint64* out = new uint64[len];
    
    uint64 pack = 0ULL;

    for(int i = 0; i < len + 15; i++){
        int tmp = 0;
        if(i < len)
            tmp = packCh(rec.seq.at(cur_headPos + i), false);
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
