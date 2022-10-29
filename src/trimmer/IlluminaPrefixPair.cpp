#include "IlluminaPrefixPair.h"

using namespace rabbit::trim;

IlluminaPrefixPair::IlluminaPrefixPair(std::string prefix1_, std::string prefix2_, rabbit::Logger& logger, int phred_, int minPrefix_, int seedMaxMiss_, int minPalindromeLikelihood_, bool palindromeKeepBoth_) : logger(logger){
    phred = phred_;
    minPrefix = minPrefix_;
    seedMaxMiss = seedMaxMiss_;
    seedMax = seedMaxMiss * 2;
    minPalindromeLikelihood = minPalindromeLikelihood_;
    palindromeKeepBoth = palindromeKeepBoth_;
    logger.infoln("Using PrefixPair: '" + prefix1_ + "' and '" + prefix2_ + "'");
    int prefix1Len = prefix1_.length();
    int prefix2Len = prefix2_.length();
    int minLength = prefix1Len < prefix2Len ? prefix1Len : prefix2Len;
    prefix1 = prefix1_.substr(0, minLength) ;
    prefix2 = prefix2_.substr(0, minLength) ;
    prefixLen = minLength;
    logger.infoln("Using PrefixPair: '" + prefix1 + "' and '" + prefix2 + "'"); // debug
    
}

int IlluminaPrefixPair::packCh(char ch, bool reverse){
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
    return 0;
}

//TODO 将Prefix 和 rec.seq 合并后再求packs, 和当前的方式比较看看哪一个效率更高
uint64* IlluminaPrefixPair::packSeqInternal(Reference& rec, bool reverse){
    int seqLen = rec.length;
    int len = prefixLen + seqLen;
    assert(len > 15);
    int cur_headPos = rec.headPos;
    uint64* out; // 用长度为16的窗口生产kmer
    out = new uint64[len - 15];
    if(!reverse){
        uint64 pack = 0;
        for(int i = 0; i < prefixLen; i++){
            int tmp = pachCh(prefix1.at(i), false);
            pack = (pack << 4) | tmp;
            if(i >= 15) out[i - 15] = pack; // 类似于长度为16的kmer
        }
        for(int i = 0; i < seqLen; i++) {
            int tmp = pachCh(rec.seq.at(i + cur_headPos), false);
            pack = (pack << 4) | tmp;
            if(i + prefixLen >= 15) out[i + prefixLen - 15] = pack; // 类似于长度为16的kmer
        }
    }
    else{
        uint64 pack = 0;
        for(int i = 0; i < prefixLen; i++){
            int tmp = pachCh(prefix2.at(i), true);
            pack = (pack >> 4) | (tmp << 60);
            if(i >= 15) out[i - 15] = pack;
        }
        for(int i = 0; i < seqLen; i++) {
            int tmp = pachCh(rec.seq.at(i + cur_headPos), true);
            pack = (pack >> 4) | (tmp << 60); // 当rec是reverse read时，产生其反向互补链的pack
            if(i + prefixLen >= 15) out[i + prefixLen - 15] = pack; // 类似于长度为16的kmer
        }
    }
    return out;
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
    if(rec1Len <= 15 || rec2Len <= 15) return std::INT_MAX; // TODO ? why

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
    return std::INT_MAX;
}