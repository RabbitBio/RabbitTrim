// #include <iomanip>
#include "trimmer/IlluminaMediumClippingSeq.h"

using namespace rabbit::trim;


IlluminaMediumClippingSeq::IlluminaMediumClippingSeq(rabbit::Logger& logger_, int phred_, std::string seq_, int seedMaxMiss_, int minSequenceLikelihood_, int minSequenceOverlap_) :logger(logger_){
    logger.infoln("Using Medium Clipping Sequence: '" + seq_ + "'");
    phred = phred_;
    seq = seq_;
    seqLen = seq.length();
    seedMaxMiss = seedMaxMiss_;
    seedMax = seedMaxMiss * 2;
    minSequenceLikelihood = minSequenceLikelihood_;
    minSequenceOverlap = minSequenceOverlap_;
    // packSeqInternal 
    pack = new uint64[seqLen - 15]; // seqLen - 16 + 1
    uint64 pack_ = 0ULL;
    for(int i = 0; i < seqLen; i++) {
        int tmp =  packCh(seq.at(i));
        pack_ = (pack_ << 4) | tmp;
        if(i >= 15) pack[i - 15] = pack_;
    }
    

}
IlluminaMediumClippingSeq::~IlluminaMediumClippingSeq(){
    delete [] pack;
}

int IlluminaMediumClippingSeq::readsSeqCompare(Reference& rec){
    std::set<int> offsetSet;
    // packSeqExternal(rec)
    uint64* packRec = packSeqExternal(rec);
    uint64* packClip = pack; 

    int packRecLen = rec.length;
    int packRecMax = packRecLen - minSequenceOverlap;
    int packClipMax = seqLen - 15; // pack.length Note：pack is produced by 'packSeqInternal'
    
    for(int i = 0; i < packRecMax; i++){
        uint64 comboMask = calcSingleMask(packRecLen - i);
        for(int j = 0; j < packClipMax; j++){
            int diff = __builtin_popcountll((packRec[i] ^ pack[j]) & comboMask);
            if(diff <= seedMax){
                int offset = i - j;
                offsetSet.emplace(offset);
            }
        }
    }

    // Iterate through offsetSet from smallest to largest
    for(auto iter = offsetSet.begin(); iter != offsetSet.end(); iter++){
        int offset = *iter;
        int recCompLength = offset > 0 ? rec.length - offset : rec.length;
        int clipCompLength = offset < 0 ? seqLen + offset : seqLen;
        int compLength = recCompLength < clipCompLength ? recCompLength : clipCompLength;
        
        assert(compLength > minSequenceOverlap);
        float seqLikelihood = calculateDifferenceQuality(rec, compLength, offset);
        if(seqLikelihood >= minSequenceLikelihood) return offset;
    }
    // return std::INT_MAX;
    return 1 << 30;
}
int IlluminaMediumClippingSeq::readsSeqCompare(neoReference& rec){
    std::set<int> offsetSet;
    // packSeqExternal(rec)
    uint64* packRec = packSeqExternal(rec);
    uint64* packClip = pack; 

    int packRecLen = rec.lseq;
    int packRecMax = packRecLen - minSequenceOverlap;
    int packClipMax = seqLen - 15; // pack.length Note：pack is produced by 'packSeqInternal'

    // std::cout << "seedMax: " << seedMax << std::endl;
    
    for(int i = 0; i < packRecMax; i++){
        uint64 comboMask = calcSingleMask(packRecLen - i);
        for(int j = 0; j < packClipMax; j++){
            int diff = __builtin_popcountll((packRec[i] ^ pack[j]) & comboMask);
            if(diff <= seedMax){
                int offset = i - j;
                offsetSet.emplace(offset);
                // std::cout << "add offset : " << offset << std::endl;
            }
            // std::cout << "i: " << i << " j: " << j << " ";
            // std::cout << "comboMask: " << std::setbase(16) << comboMask << " ";
            // std::cout << "packRec: " << packRec[i] << " ";
            // std::cout << "packClip: " << packClip[j] << " ";
            // std::cout << "diff: " << std::setbase(10) << diff << std::endl;
        }
    }

    // Iterate through offsetSet from smallest to largest
    for(auto iter = offsetSet.begin(); iter != offsetSet.end(); iter++){
        int offset = *iter;
        int recCompLength = offset > 0 ? rec.lseq - offset : rec.lseq;
        int clipCompLength = offset < 0 ? seqLen + offset : seqLen;
        int compLength = recCompLength < clipCompLength ? recCompLength : clipCompLength;
        
        assert(compLength > minSequenceOverlap);
        float seqLikelihood = calculateDifferenceQuality(rec, compLength, offset);
        if(seqLikelihood >= minSequenceLikelihood) return offset;
    }
    // return std::INT_MAX;
    return 1 << 30;
}


int IlluminaMediumClippingSeq::packCh(char ch){
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

uint64* IlluminaMediumClippingSeq::packSeqExternal(Reference& rec){ 
    int len = rec.length;
    int cur_headPos = rec.headPos;
    uint64* out = new uint64[len];
    
    uint64 pack = 0ULL;

    for(int i = 0; i < len + 15; i++){
        uint64 tmp = 0;
        if(i < len)
            tmp = packCh(rec.seq.at(cur_headPos + i));
        pack = (pack << 4) | tmp;
        if(i >= 15) out[i - 15] = pack;
    }
    return out;
}
uint64* IlluminaMediumClippingSeq::packSeqExternal(neoReference& rec){ 
    // TODO Vectorize
    int len = rec.lseq;
    char* rec_seq = (char*) (rec.base + rec.pseq);
    uint64* out = new uint64[len];
    uint64 pack = 0ULL;

    for(int i = 0; i < len + 15; i++){
        uint64 tmp = 0;
        if(i < len)
            tmp = packCh(rec_seq[i]);
        pack = (pack << 4) | tmp;
        if(i >= 15) out[i - 15] = pack;
    }
    return out;
}

uint64 IlluminaMediumClippingSeq::calcSingleMask(int length){
    uint64 mask = MASK_VAL;
    if(length < 16) {
        mask <<= (16 - length) * 4;
    }
    return mask;
}

float IlluminaMediumClippingSeq::calculateMaximumRange(float* vals, int valsLen){
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

float IlluminaMediumClippingSeq::calculateDifferenceQuality(Reference& rec, int overlap, int recOffset){
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

float IlluminaMediumClippingSeq::calculateDifferenceQuality(neoReference& rec, int overlap, int recOffset){
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

