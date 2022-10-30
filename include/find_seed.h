#include "func_set.h"
#include "array_set.h"
#include <cstdint>
#include <immintrin.h>
using namespace std;

// avx512 左移8bit
__m512i idx =  _mm512_set_epi32(11,10,9,8,7,6,5,4,3,2,1,0,0,0,0,0);
inline __m512i _mm512_fake_slli_si512_1(__m512i a, __m512i  zero_, __m512i idx){
    __m512i sl = _mm512_mask_permutexvar_epi32(zero_, 0xfff0,  idx, a);
    return _mm512_alignr_epi8(a,sl,(16-1));
    //return _mm512_alignr_epi8(a,sl,16-1);
}

void find_seed(char* seq,int num_, const int16_t** seed_pos_1, const int16_t** seed_pos_2, int p1, int p2,
                char index_10,char index_11,char index_12,char index_20,char index_21,char index_22){
    // num_ = len(seq) / 62;
    // if (len(seq) % 62 > 2) num_++

    __m512i m_index_10 = _mm512_set1_epi8 (index_12);
    __m512i m_index_11 = _mm512_set1_epi8 (index_11);
    __m512i m_index_12 = _mm512_set1_epi8 (index_10);
    __m512i m_index_20 = _mm512_set1_epi8 (index_22);
    __m512i m_index_21 = _mm512_set1_epi8 (index_21);
    __m512i m_index_22 = _mm512_set1_epi8 (index_20);

    __m512i  zero_ = _mm512_setzero_si512();
    for(int i = 0; i < num_; i++){
        __m512i load1 =  _mm512_loadu_si512(seq + ((i << 6) - (i << 1))); // seq + 62i
        __m512i res_10 = _mm512_xor_si512 (load1, m_index_10);
        __m512i res_20 = _mm512_xor_si512 (load1, m_index_20);
        // 左移8bit
        load1 = _mm512_fake_slli_si512_1(load1, zero_, idx);
        __m512i res_11 = _mm512_xor_si512 (load1, m_index_11);
        res_11 = _mm512_or_si512 (res_10, res_11);
        __m512i res_21 = _mm512_xor_si512 (load1, m_index_21);
        res_21 = _mm512_or_si512 (res_20, res_21);

        load1 = _mm512_fake_slli_si512_1(load1, zero_, idx);
        __m512i res_12 = _mm512_xor_si512 (load1, m_index_12);
        res_12 = _mm512_or_si512 (res_11, res_12);
        __mmask64 res_1 = _mm512_cmp_epu8_mask(res_12,zero_,0);
        unsigned long long int64_res_1 = res_1;
        __m512i res_22 = _mm512_xor_si512 (load1, m_index_22);
        res_22 = _mm512_or_si512 (res_21, res_22);
        //__mmask64 res = _mm512_mask_cmp_epu8_mask(0xffffffffffffffff,res_12,zero_,0);
        __mmask64 res_2 = _mm512_cmp_epu8_mask(res_22,zero_,0);
        unsigned long long int64_res_2 = res_2;
        
        // const int16_t * seed_arr_1;
        // const int16_t * seed_arr_2;
        for(int j = 0; j < 8; j++){
            int t_1  = int64_res_1 & 255;
            int t_2  = int64_res_2 & 255;
            int64_res_1 = int64_res_1 >> 8;
            int64_res_2 = int64_res_2 >> 8;
            // 确定种子的位置
            seed_pos_1[p1++] = seed_table[i][j][t_1];
            seed_pos_2[p2++] = seed_table[i][j][t_2];
            // seed_arr_1 = seed_table[i][j][t_1];
            // seed_arr_2 = seed_table[i][j][t_2]; // -3
        }
    }
}
