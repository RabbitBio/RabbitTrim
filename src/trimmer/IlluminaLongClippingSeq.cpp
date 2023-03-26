#include "trimmer/IlluminaLongClippingSeq.h"
#include <iomanip>
#if defined __SSE2__ && defined __AVX__ && defined __AVX2__ && defined TRIM_USE_VEC
#include <immintrin.h>
#endif

using namespace rabbit::trim;

IlluminaLongClippingSeq::IlluminaLongClippingSeq(rabbit::Logger& logger_, int phred_, std::string seq_, int seedMaxMiss_, int minSequenceLikelihood_, int minSequenceOverlap_, int consumerNum_) : logger(logger_){
    logger.infoln("Using Long Clipping Sequence: '" + seq_ + "'");
    phred = phred_;
    seq = seq_;
    seqLen = seq.length();
    seedMaxMiss = seedMaxMiss_;
    seedMax = seedMaxMiss * 2;
    minSequenceLikelihood = minSequenceLikelihood_;
    minSequenceOverlap = minSequenceOverlap_;
    consumerNum = consumerNum_;
    
    // packSeqInternal(seq)
    fullPack = new uint64[seqLen - 15]; // seqLen - 16 + 1
    uint64 pack_ = 0ULL;
    for(int i = 0; i < seqLen; i++){
        uint64 tmp =  (1 << ((seq[i] >> 1) & 7)) & 15;
        pack_ = (pack_ << 4) | tmp;
        if(i >= 15) fullPack[i - 15] = pack_;
    }
    pack = new uint64[(seqLen - 15 + INTERLEAVE - 1) / INTERLEAVE];
    for(int i = 0; i < seqLen - 15; i += INTERLEAVE){
        pack[i / INTERLEAVE] = fullPack[i];
    }

    recPacks = new uint64[(rabbit::trim::MAX_READ_LENGTH + 15) * consumerNum_];
    likelihoodTotal = new float[(rabbit::trim::MAX_READ_LENGTH) * consumerNum_];
#if defined __SSE2__ && defined __AVX__ && defined __AVX2__ && defined TRIM_USE_VEC
    logger_.errorln("Using vectorized in Long Clipping ...");
    
    seq_str = new char[seqLen + 31];
    for(int i = 0; i < seqLen; i++)
    {
      seq_str[i] = seq[i];
    }

    phred_arr = new char[16];
    all_N = new char[16];
    for(int i = 0; i < 16; i++)
    {
      phred_arr[i] = phred;
      all_N[i] = 'N';
    }

    awards = new float[8];
    divide_arr = new float[8];
    for(int i = 0; i < 8; i++)
    {
      awards[i] = LOG10_4;
      divide_arr[i] = -10.0f;
    }
    
#endif

}
IlluminaLongClippingSeq::~IlluminaLongClippingSeq(){
    delete [] fullPack;
    delete [] pack;
    delete [] recPacks;
    delete [] likelihoodTotal;
#if defined __SSE2__ && defined __AVX__ && defined __AVX2__ && TRIM_USE_VEC
    delete [] seq_str;
    delete [] awards;
    delete [] phred_arr;
    delete [] divide_arr;
    delete [] all_N;
#endif
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
            if(diff < seedMax){
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
    return 1 << 30;
}

int IlluminaLongClippingSeq::readsSeqCompare(neoReference& rec){
    std::set<int> offsetSet;
    uint64* packRec = packSeqExternal(rec);
    uint64* packClip = pack;

    // 根据minSequenceOverlap计算最大比较次数
    int packRecMax = rec.lseq - minSequenceOverlap;
    int packClipMax = (seqLen - 15 + INTERLEAVE - 1) / INTERLEAVE; // pack.length
    
    for(int i = 0; i < packRecMax; i++){
        uint64 comboMask = calcSingleMask(rec.lseq - i);
        for(int j = 0; j < packClipMax; j++){
            int diff = __builtin_popcountll((packRec[i] ^ packClip[j]) & comboMask);
						if(diff <= seedMax){
							offsetSet.emplace(i - j * INTERLEAVE);
						}
				}
		}

		for(auto iter = offsetSet.begin(); iter != offsetSet.end(); iter++){
			int offset = *iter;
			int recCompLength = offset > 0 ? rec.lseq - offset : rec.lseq;
			int clipCompLength = offset < 0 ? seqLen + offset : seqLen;
			int compLength = recCompLength < clipCompLength ? recCompLength : clipCompLength;

				assert(compLength > minSequenceOverlap);
			float seqLikelihood = calculateDifferenceQuality(rec, compLength, offset);
			if(seqLikelihood >= minSequenceLikelihood) return offset;
		}
		return 1 << 30;
}

int IlluminaLongClippingSeq::readsSeqCompare(neoReference& rec, int threadId){
    // std::set<int> offsetSet;
    uint64* packRec = packSeqExternal(rec, threadId);
    uint64* packClip = pack;

    // 根据minSequenceOverlap计算最大比较次数
    int packRecMax = rec.lseq - minSequenceOverlap;
    int packClipMax = (seqLen - 15 + INTERLEAVE - 1) / INTERLEAVE; // pack.length
    
    // for(int i = 0; i < packRecMax; i++){
    //     uint64 comboMask = calcSingleMask(rec.lseq - i);
    //     for(int j = 0; j < packClipMax; j++){
    //         int diff = __builtin_popcountll((packRec[i] ^ packClip[j]) & comboMask);
		// 				if(diff <= seedMax){
		// 					offsetSet.emplace(i - j * INTERLEAVE);
		// 				}

		// 		}
		// }

		// for(auto iter = offsetSet.begin(); iter != offsetSet.end(); iter++){
		// 	int offset = *iter;
		// 	int recCompLength = offset > 0 ? rec.lseq - offset : rec.lseq;
		// 	int clipCompLength = offset < 0 ? seqLen + offset : seqLen;
		// 	int compLength = recCompLength < clipCompLength ? recCompLength : clipCompLength;

		// 		assert(compLength > minSequenceOverlap);
		// 	float seqLikelihood = calculateDifferenceQuality(rec, compLength, offset);
		// 	if(seqLikelihood >= minSequenceLikelihood) return offset;
		// }
    int packRecLen = rec.lseq;
    int minOffset = 0 - (packClipMax - 1) * INTERLEAVE;
    int maxOffset = packRecMax - 1;
    int offsetNums = maxOffset - minOffset + 1;
    int cntLen = (offsetNums + 63) / 64;
    uint64* cnt = new uint64[cntLen];
    for(int i = 0; i < cntLen; i++)
    {
      cnt[i] = 0ULL;
    }
    for(int i = 0; i < packRecMax; i++){
        uint64 comboMask = calcSingleMask(packRecLen - i);
        for(int j = 0; j < packClipMax / 4 * 4; j+=4){
            int diff1 = __builtin_popcountll((packRec[i] ^ pack[j + 0]) & comboMask);
            int diff2 = __builtin_popcountll((packRec[i] ^ pack[j + 1]) & comboMask);
            int diff3 = __builtin_popcountll((packRec[i] ^ pack[j + 2]) & comboMask);
            int diff4 = __builtin_popcountll((packRec[i] ^ pack[j + 3]) & comboMask);
            if(diff1 <= seedMax){
                int offset = i - j * INTERLEAVE;
                cnt[(offset - minOffset) >> 6] |= (1ULL << (63 - ((offset - minOffset) % 64)));
            }
            if(diff2 <= seedMax){
                int offset = i - (j + 1) * INTERLEAVE;
                cnt[(offset - minOffset) >> 6] |= (1ULL << (63 - ((offset - minOffset) % 64)));
            }
            if(diff3 <= seedMax){
                int offset = i - (j + 2) * INTERLEAVE;
                cnt[(offset - minOffset) >> 6] |= (1ULL << (63 - ((offset - minOffset) % 64)));
            }
            if(diff4 <= seedMax){
                int offset = i - (j + 3) * INTERLEAVE;
                cnt[(offset - minOffset) >> 6] |= (1ULL << (63 - ((offset - minOffset) % 64)));
            }
        }
    }

    for(int i = 0; i < packRecMax; i++){
        uint64 comboMask = calcSingleMask(packRecLen - i);
        for(int j = packClipMax / 4 * 4; j < packClipMax; j++){
            int diff = __builtin_popcountll((packRec[i] ^ pack[j]) & comboMask);
            if(diff <= seedMax){
                int offset = i - j * INTERLEAVE;
                cnt[(offset - minOffset) >> 6] |= (1ULL << (63 - ((offset - minOffset) % 64)));
            }
        }
    }

    for(int i = 0; i < cntLen; i++)
    {
      for(int p = 0; p < 64; p++)
      {
        if((cnt[i] >> (63 - p)) & 1)
        {
          int offset = minOffset + i * 64 + p;
          int recCompLength = offset > 0 ? rec.lseq - offset : rec.lseq;
          int clipCompLength = offset < 0 ? seqLen + offset : seqLen;
          int compLength = recCompLength < clipCompLength ? recCompLength : clipCompLength;

          assert(compLength > minSequenceOverlap);
          float seqLikelihood = calculateDifferenceQuality(rec, compLength, offset, threadId);
          if(seqLikelihood >= minSequenceLikelihood) return offset;
        }
      }
    }
		return 1 << 30;
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
		uint64 tmp = 0;
		if(i < len)
			tmp = packCh(rec.seq.at(cur_headPos + i));
		pack = (pack << 4) | tmp;
		if(i >= 15) out[i - 15] = pack;
	}
	return out;
}
uint64* IlluminaLongClippingSeq::packSeqExternal(neoReference& rec){ 
	int len = rec.lseq;
	char* rec_seq = (char*)(rec.base + rec.pseq);
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

uint64* IlluminaLongClippingSeq::packSeqExternal(neoReference& rec, int threadId){ 
	int len = rec.lseq;
	char* rec_seq = (char*)(rec.base + rec.pseq);
	uint64* out = recPacks + threadId * (rabbit::trim::MAX_READ_LENGTH + 15); 
	uint64 pack = 0ULL;
  
  for(int i = 0; i < len / 4 * 4; i+=4)
  {
    uint64 tmp0 = (1 << ((rec_seq[i + 0] >> 1) & 7) ) & 15;
    uint64 tmp1 = (1 << ((rec_seq[i + 1] >> 1) & 7) ) & 15;
    uint64 tmp2 = (1 << ((rec_seq[i + 2] >> 1) & 7) ) & 15;
    uint64 tmp3 = (1 << ((rec_seq[i + 3] >> 1) & 7) ) & 15;
		pack = (pack << 4) | tmp0;
		out[i + 0] = pack;
		pack = (pack << 4) | tmp1;
		out[i + 1] = pack;
		pack = (pack << 4) | tmp2;
		out[i + 2] = pack;
		pack = (pack << 4) | tmp3;
		out[i + 3] = pack;
  }
  for(int i = len / 4 * 4; i < len; i++)
  {
    uint64 tmp = (1 << ((rec_seq[i] >> 1) & 7) ) & 15;
		pack = (pack << 4) | tmp;
		out[i] = pack;
  }

	for(int i = len; i < len + 15; i+=3){
		uint64 pack1 = pack << 4;
		uint64 pack2 = pack << 8;
		uint64 pack3 = pack << 12;
    pack = pack3;
		out[i + 0] = pack1;
		out[i + 1] = pack2;
		out[i + 2] = pack3;
    
	}
	return out + 15;
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
float IlluminaLongClippingSeq::calculateDifferenceQuality(neoReference& rec, int overlap, int recOffset){
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
        char ch2 = seq[clipPos];
        
        int qual_val = rec_qual[recPos] - phred;
        float s = ((ch1 >> 1) & 3) == ((ch2 >> 1) & 3) ? LOG10_4 : -qual_val / 10.0f;
        likelihood[i] = ((ch1 == 'N' || ch2 == 'N') ? 0 : s);
        recPos++;
				clipPos++;
		}
		float l = calculateMaximumRange(likelihood, overlap);
		return l;
}

float IlluminaLongClippingSeq::calculateDifferenceQuality(neoReference& rec, int overlap, int recOffset, int threadId){
  // getQualityAsInteger
  int len = rec.lseq;
  char* rec_seq = (char*)(rec.base + rec.pseq);
  char* rec_qual = (char*)(rec.base + rec.pqual);

  // Location to start comparison
  int recPos = recOffset > 0 ? recOffset : 0;
  int clipPos = recOffset < 0 ? -recOffset : 0;
  float* likelihood = likelihoodTotal + rabbit::trim::MAX_READ_LENGTH * threadId;
#if defined __SSE2__ && defined __AVX__ && defined __AVX2__ && defined TRIM_USE_VEC
  // 16B extra space required for rear
  char* tmp_rec_qual = rec_qual + recPos;
  char* tmp_rec_seq = rec_seq + recPos;
  char* tmp_seq_str = seq_str + clipPos;
  int nums = (overlap + 15) / 16;
  
  for(int i = 0; i < nums; i++)
  {
    __m128i vphred  = _mm_loadu_si128((__m128i_u*)phred_arr);
    __m128i vqual  = _mm_loadu_si128((__m128i_u*)(tmp_rec_qual + i * 16));
    __m128i vscore  = _mm_sub_epi8(vqual, vphred);
    __m256i vscore_1 = _mm256_cvtepi8_epi32(vscore);
    __m128i vscore_sr8 = _mm_srli_si128(vscore, 0x8);
    __m256i vscore_2 = _mm256_cvtepi8_epi32(vscore_sr8);
    __m256 vscore_1_ps = _mm256_cvtepi32_ps(vscore_1);
    __m256 vscore_2_ps = _mm256_cvtepi32_ps(vscore_2);
    __m256 vdivide = _mm256_loadu_ps(divide_arr);
    __m256 vscore_0 = _mm256_setzero_ps();
    vscore_1_ps = _mm256_div_ps(vscore_1_ps, vdivide);
    vscore_2_ps = _mm256_div_ps(vscore_2_ps, vdivide);

    // vscore_1_ps = _mm256_ceil_ps(vscore_1_ps);
    // vscore_2_ps = _mm256_ceil_ps(vscore_2_ps);

    __m128i vrec  = _mm_loadu_si128((__m128i_u*)(tmp_rec_seq + i * 16));

    __m128i vallN = _mm_loadu_si128((__m128i_u*)(all_N));
    __m128i vres3  = _mm_cmpeq_epi8(vrec, vallN);
    __m128i vres3_sr8 = _mm_srli_si128(vres3, 0x8);
    __m256i vlow3 = _mm256_cvtepi8_epi32(vres3);
    __m256i vhigh3 = _mm256_cvtepi8_epi32(vres3_sr8);
    __m256 vlow3_ps = _mm256_cvtepi32_ps(vlow3);
    __m256 vhigh3_ps = _mm256_cvtepi32_ps(vhigh3);
    vscore_1_ps = _mm256_blendv_ps(vscore_1_ps, vscore_0, vlow3_ps);
    vscore_2_ps = _mm256_blendv_ps(vscore_2_ps, vscore_0, vhigh3_ps);
    __m128i vclip = _mm_loadu_si128((__m128i_u*)(tmp_seq_str + i * 16));
    __m128i vres1  = _mm_cmpeq_epi8(vrec, vclip);
    __m128i vres1_sr8 = _mm_srli_si128(vres1, 0x8);
    __m256i vlow1 = _mm256_cvtepi8_epi32(vres1);
    __m256i vhigh1 = _mm256_cvtepi8_epi32(vres1_sr8);
    __m256 vlow1_ps = _mm256_cvtepi32_ps(vlow1);
    __m256 vhigh1_ps = _mm256_cvtepi32_ps(vhigh1);
    __m256 vawards = _mm256_loadu_ps(awards);
    vscore_1_ps = _mm256_blendv_ps(vscore_1_ps, vawards, vlow1_ps);
    vscore_2_ps = _mm256_blendv_ps(vscore_2_ps, vawards, vhigh1_ps);
    _mm256_storeu_ps(likelihood + i * 16,      vscore_1_ps);
    _mm256_storeu_ps(likelihood + i * 16 + 8,  vscore_2_ps);
  }

#else
  for(int i = 0; i < overlap; i++){
    char ch1 = rec_seq[recPos + i];
    char ch2 = seq[clipPos + i];

    float penalty = -1 * ((bool)((1 << ((ch1 >> 1) & 7)) & 15)) * (rec_qual[recPos + i] - phred) / 10.0f;
    float s = (((ch1 >> 1) & 3) == ((ch2 >> 1) & 3)) * LOG10_4 + (((ch1 >> 1) & 3) != ((ch2 >> 1) & 3)) * penalty;
    likelihood[i] = s;
  }
#endif
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
