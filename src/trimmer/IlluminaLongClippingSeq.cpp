#include "IlluminaLongClippingSeq.h"

using namespace rabbit;


IlluminaLongClippingSeq::IlluminaLongClippingSeq(rabbit::Logger& logger_, int phred_, std::string seq_, int seedMaxMiss_, int minSequenceLikelihood_, int minSequenceOverlap_) : logger(logger_){
    logger.infoln("Using Long Clipping Sequence: '" + seq_ + "'");
    phred = phred_;
    seq = seq_;
    seqLen = seq.length();
    seedMaxMiss = seedMaxMiss_;
    seedMax = seedMaxMiss * 2;
    minSequenceLikelihood = minSequenceLikelihood_;
    minSequenceOverlap = minSequenceOverlap_;
    
    // packSeqInternal(seq)
    fullPack = new uint64[seqLen - 15]; // seqLen - 16 + 1
    uint64 pack_ = 0ULL;
    for(int i = 0; i < seqLen; i++){
        char ch = seq.at(i);
        int tmp = packCh(ch);
        pack_ = (pack_ << 4) | tmp;
        if(i >= 15) fullPack[i - 15] = pack_;
    }
    pack = new uint64[(seqLen - 15 + INTERLEAVE - 1) / INTERLEAVE];
    for(int i = 0; i < seqLen - 15; i += INTERLEAVE){
        pack[i / INTERLEAVE] = fullPack[i];
    }

}
IlluminaLongClippingSeq::~IlluminaLongClippingSeq(){
    delete [] fullPack;
    delete [] pack;
}

int IlluminaLongClippingSeq::readsSeqCompare(Reference& rec){
    std::set<int> offsetSet;
    uint64* packRec = packSeqExternal(rec);
    uint64* packClip = pack;

    // 根据minSequenceOverlap计算最大比较次数
    int packRecMax = rec.length - minSequenceOverlap;
    int packClipMax = (seqLen - 15 + INTERLEAVE - 1) / INTERLEAVE; // pack.length
    
    for(int i = 0; i < packRecMax; i++){
        uint64 comboMask = calcSingleMask(rec.length - i);
        for(int j = 0; j < packClipMax; j++){
            int diff = __builtin_popcountll((packRec[i] ^ packClip[j]) & comboMask);
            if(diff <= seedMax){
                offsetSet.emplace(i - j * INTERLEAVE);
            }
        }
    }
    
    for(auto iter = offsetSet.begin(); iter != offsetSet.end(); iter++){
        int offset = *iter;
        int recCompLength = offset > 0 ? rec.length - offset : rec.length;
        int clipCompLength = offset < 0 ? seqLen + offset : seqLen;
        int compLength = recCompLength < clipCompLength ? recCompLength : clipCompLength;
        
        assert(compLength > minSequenceOverlap);
        float seqLikelihood = calculateDifferenceQuality(rec, compLength, offset);
        if(seqLikelihood >= minSequenceLikelihood) return offset;
    }
    return std::INT_MAX;
}

int IlluminaLongClippingSeq::packCh(char ch){
    switch (ch)
    {
    case 'A':
        return BASE_A;
    case 'T':
        return BASE_T;
    case 'C':
        return BASE_C;
    case 'G':
        return BASE_G;
    }
    return 0;
}


uint64* IlluminaLongClippingSeq::packSeqExternal(Reference& rec){ 
    int len = rec.length;
    int cur_headPos = rec.headPos;
    uint64* out = new uint64[len];
    
    uint64 pack = 0ULL;

    for(int i = 0; i < len + 15; i++){
        int tmp = 0;
        if(i < seqLen)
            tmp = packCh(rec.seq.at(cur_headPos + i), false);
        pack = (pack << 4) | tmp;
        if(i >= 15) out[i - 15] = pack;
    }
    return out;
}


uint64 IlluminaLongClippingSeq::calcSingleMask(int length){
    uint64 mask = MASK_VAL;
    if(length < 16) {
        mask <<= (16 - length) * 4;
    }
    return mask;
}

float IlluminaLongClippingSeq::calculateDifferenceQuality(Reference& rec, int overlap, int recOffset){
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


float IlluminaLongClippingSeq::calculateMaximumRange(float* vals, int valsLen){
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