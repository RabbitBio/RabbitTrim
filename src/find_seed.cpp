#include <iostream>
#include <string.h>
#include <immintrin.h>
#include "stdio.h"
#include "stdlib.h"
#include "array_set.h"
#include <cstdint>
#include "util.h"


using namespace std;
// avx512 左移8bit
__m512i idx =  _mm512_set_epi32(11,10,9,8,7,6,5,4,3,2,1,0,0,0,0,0);
inline __m512i _mm512_fake_slli_si512_1(__m512i a, __m512i  zero_, __m512i idx){
    __m512i sl = _mm512_mask_permutexvar_epi32(zero_, 0xfff0,  idx, a);
    return _mm512_alignr_epi8(a,sl,(16-1));
}

inline void print_avs512_seq(__m512i& v){
    char seq_1[65];
    seq_1[64] = '\0';
    _mm512_storeu_si512((__m512i*)seq_1, v);
    for(int i = 0; i < 64; i++){
        printf("%c", seq_1[i]); 
    }
    printf("\n");
    fflush(stdout);
}


void find_seed(char* seq,int num, const int16_t** seed_pos_1, const int16_t** seed_pos_2, int p1, int p2,
                char index_10,char index_11,char index_12,char index_20,char index_21,char index_22){

    // seq_len  = 62 * n + 2
    //int n = len / 62;
    //if(n % 62 > 2) n++;

    // index1 = "CTG"
    //__m512i m_index_10 = _mm512_set1_epi8 (index_12);
    //__m512i m_index_11 = _mm512_set1_epi8 (index_11);
    //__m512i m_index_12 = _mm512_set1_epi8 (index_10);
    //__m512i m_index_20 = _mm512_set1_epi8 (index_22);
    //__m512i m_index_21 = _mm512_set1_epi8 (index_21);
    //__m512i m_index_22 = _mm512_set1_epi8 (index_20);
    __m512i m_index_10 = _mm512_set1_epi8 (index_10);
    __m512i m_index_11 = _mm512_set1_epi8 (index_11);
    __m512i m_index_12 = _mm512_set1_epi8 (index_12);
    __m512i m_index_20 = _mm512_set1_epi8 (index_20);
    __m512i m_index_21 = _mm512_set1_epi8 (index_21);
    __m512i m_index_22 = _mm512_set1_epi8 (index_22);


    __m512i  zero_ = _mm512_setzero_si512();
#pragma unroll(4)
    for(int i = 0; i < num; i++){
        //__m512i load1 =  _mm512_loadu_si512(seq + ((i << 6) - (i << 1))); // seq + 62i
        __m512i load1 =  _mm512_loadu_si512(seq + (i << 6)); // seq + 64i
        __m512i res_10 = _mm512_xor_si512 (load1, m_index_10);
        __m512i res_20 = _mm512_xor_si512 (load1, m_index_20);
        // 左移8bit
        //load1 = _mm512_fake_slli_si512_1(load1, zero_, idx);
        load1 = _mm512_loadu_si512(seq + ((i << 6) + (1)));
        __m512i res_11 = _mm512_xor_si512 (load1, m_index_11);
        res_11 = _mm512_or_si512 (res_10, res_11);
        __m512i res_21 = _mm512_xor_si512 (load1, m_index_21);
        res_21 = _mm512_or_si512 (res_20, res_21);

        //load1 = _mm512_fake_slli_si512_1(load1, zero_, idx);
        load1 = _mm512_loadu_si512(seq + ((i << 6) + (2)));
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

int main(){
    const char* file_name = "./testing_dataset/isReads.read1.fq";
    FILE* fp = NULL;
    char * seq = new char[622];
    fp = fopen(file_name,"r");
    if(!fp){
        cout << "error "<<endl;
    }
    double start,finish;
    double total_time = 0;

    for (int i = 0; i < (1 << 23); i++){
        memset(seq,0,622);
        fgets(seq,622,fp); // name
        fgets(seq,622,fp); // seq
        int len = strlen(seq) -1;
        const int16_t ** seed_pos_1 = new const int16_t*[72]; // 9*8
        const int16_t ** seed_pos_2 = new const int16_t*[72]; // 9*8
        for(int i = 0; i< 72;i++){
            seed_pos_1[i] = NULL;
            seed_pos_2[i] = NULL;
        }
        
        int num = len   / 64;
        if((len % 64) > 0){
            num++;
        }
        //int num = (len -2) / 62;
        //if((len % 62) > 2){
        //    num++;
        //}
        start = getTime();
        find_seed(seq,num,seed_pos_1, seed_pos_2, 0, 0, 'C','T','G','T','C','T');
        // register char *poffset  =seq;
	    // register char *indexloc = poffset;
	    // while( true ) {
		//     indexloc = strstr( indexloc, "CTG" ); // strstr(char*a, char*b)是在a中找到第一次出现b的位置 返回值是char* 找不到时返回NULL
		//     if( indexloc == NULL )
		// 	    break;
		//     seed[( indexloc - poffset )] = 1; // 
		//     indexloc ++;
	    // }
        finish = getTime();
        // cout<<seq;
        //for(int i = 0; i< 72;i++){
        //    if(seed_pos_1[i]==NULL && seed_pos_2[i] == NULL)
        //        break;
        //    if(seed_pos_1[i]){
        //        for(int j = 0; j < seed_pos_1[i][0];j++){
        //            // cout << seed_pos_1[i][j+1] << " ";
        //            seed_pos_1[i][j+1];
        //        }
        //    }
        //    if(seed_pos_2[i]){
        //        for(int j = 0; j < seed_pos_2[i][0];j++){
        //            // cout << seed_pos_2[i][j+1]-3 << " ";
        //            seed_pos_2[i][j+1]-3;
        //        }
        //    }
        //}
        // cout << endl;
        delete [] seed_pos_1;
        delete [] seed_pos_2;
        fgets(seq,622,fp); // strand
        fgets(seq,622,fp); // quality
        total_time += finish -start;
    }
    fclose(fp);
    
    cout<<"Time : " << total_time << " s" <<endl;
    return 0;
}

