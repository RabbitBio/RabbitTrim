#include <vector>
using namespace std;



inline void func0( int* seed, int i, int j ){
}
inline void func1( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+7] = 1;
}
inline void func2( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+6] = 1;
}
inline void func3( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+6] = 1;
   seed[p+7] = 1;
}
inline void func4( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+5] = 1;
}
inline void func5( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+5] = 1;
   seed[p+7] = 1;
}
inline void func6( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+5] = 1;
   seed[p+6] = 1;
}
inline void func7( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+5] = 1;
   seed[p+6] = 1;
   seed[p+7] = 1;
}
inline void func8( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+4] = 1;
}
inline void func9( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+4] = 1;
   seed[p+7] = 1;
}
inline void func10( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+4] = 1;
   seed[p+6] = 1;
}
inline void func11( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+4] = 1;
   seed[p+6] = 1;
   seed[p+7] = 1;
}
inline void func12( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+4] = 1;
   seed[p+5] = 1;
}
inline void func13( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+4] = 1;
   seed[p+5] = 1;
   seed[p+7] = 1;
}
inline void func14( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+4] = 1;
   seed[p+5] = 1;
   seed[p+6] = 1;
}
inline void func15( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+4] = 1;
   seed[p+5] = 1;
   seed[p+6] = 1;
   seed[p+7] = 1;
}
inline void func16( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+3] = 1;
}
inline void func17( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+3] = 1;
   seed[p+7] = 1;
}
inline void func18( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+3] = 1;
   seed[p+6] = 1;
}
inline void func19( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+3] = 1;
   seed[p+6] = 1;
   seed[p+7] = 1;
}
inline void func20( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+3] = 1;
   seed[p+5] = 1;
}
inline void func21( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+3] = 1;
   seed[p+5] = 1;
   seed[p+7] = 1;
}
inline void func22( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+3] = 1;
   seed[p+5] = 1;
   seed[p+6] = 1;
}
inline void func23( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+3] = 1;
   seed[p+5] = 1;
   seed[p+6] = 1;
   seed[p+7] = 1;
}
inline void func24( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+3] = 1;
   seed[p+4] = 1;
}
inline void func25( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+3] = 1;
   seed[p+4] = 1;
   seed[p+7] = 1;
}
inline void func26( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+3] = 1;
   seed[p+4] = 1;
   seed[p+6] = 1;
}
inline void func27( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+3] = 1;
   seed[p+4] = 1;
   seed[p+6] = 1;
   seed[p+7] = 1;
}
inline void func28( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+3] = 1;
   seed[p+4] = 1;
   seed[p+5] = 1;
}
inline void func29( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+3] = 1;
   seed[p+4] = 1;
   seed[p+5] = 1;
   seed[p+7] = 1;
}
inline void func30( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+3] = 1;
   seed[p+4] = 1;
   seed[p+5] = 1;
   seed[p+6] = 1;
}
inline void func31( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+3] = 1;
   seed[p+4] = 1;
   seed[p+5] = 1;
   seed[p+6] = 1;
   seed[p+7] = 1;
}
inline void func32( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+2] = 1;
}
inline void func33( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+2] = 1;
   seed[p+7] = 1;
}
inline void func34( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+2] = 1;
   seed[p+6] = 1;
}
inline void func35( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+2] = 1;
   seed[p+6] = 1;
   seed[p+7] = 1;
}
inline void func36( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+2] = 1;
   seed[p+5] = 1;
}
inline void func37( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+2] = 1;
   seed[p+5] = 1;
   seed[p+7] = 1;
}
inline void func38( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+2] = 1;
   seed[p+5] = 1;
   seed[p+6] = 1;
}
inline void func39( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+2] = 1;
   seed[p+5] = 1;
   seed[p+6] = 1;
   seed[p+7] = 1;
}
inline void func40( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+2] = 1;
   seed[p+4] = 1;
}
inline void func41( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+2] = 1;
   seed[p+4] = 1;
   seed[p+7] = 1;
}
inline void func42( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+2] = 1;
   seed[p+4] = 1;
   seed[p+6] = 1;
}
inline void func43( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+2] = 1;
   seed[p+4] = 1;
   seed[p+6] = 1;
   seed[p+7] = 1;
}
inline void func44( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+2] = 1;
   seed[p+4] = 1;
   seed[p+5] = 1;
}
inline void func45( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+2] = 1;
   seed[p+4] = 1;
   seed[p+5] = 1;
   seed[p+7] = 1;
}
inline void func46( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+2] = 1;
   seed[p+4] = 1;
   seed[p+5] = 1;
   seed[p+6] = 1;
}
inline void func47( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+2] = 1;
   seed[p+4] = 1;
   seed[p+5] = 1;
   seed[p+6] = 1;
   seed[p+7] = 1;
}
inline void func48( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+2] = 1;
   seed[p+3] = 1;
}
inline void func49( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+2] = 1;
   seed[p+3] = 1;
   seed[p+7] = 1;
}
inline void func50( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+2] = 1;
   seed[p+3] = 1;
   seed[p+6] = 1;
}
inline void func51( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+2] = 1;
   seed[p+3] = 1;
   seed[p+6] = 1;
   seed[p+7] = 1;
}
inline void func52( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+2] = 1;
   seed[p+3] = 1;
   seed[p+5] = 1;
}
inline void func53( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+2] = 1;
   seed[p+3] = 1;
   seed[p+5] = 1;
   seed[p+7] = 1;
}
inline void func54( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+2] = 1;
   seed[p+3] = 1;
   seed[p+5] = 1;
   seed[p+6] = 1;
}
inline void func55( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+2] = 1;
   seed[p+3] = 1;
   seed[p+5] = 1;
   seed[p+6] = 1;
   seed[p+7] = 1;
}
inline void func56( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+2] = 1;
   seed[p+3] = 1;
   seed[p+4] = 1;
}
inline void func57( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+2] = 1;
   seed[p+3] = 1;
   seed[p+4] = 1;
   seed[p+7] = 1;
}
inline void func58( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+2] = 1;
   seed[p+3] = 1;
   seed[p+4] = 1;
   seed[p+6] = 1;
}
inline void func59( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+2] = 1;
   seed[p+3] = 1;
   seed[p+4] = 1;
   seed[p+6] = 1;
   seed[p+7] = 1;
}
inline void func60( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+2] = 1;
   seed[p+3] = 1;
   seed[p+4] = 1;
   seed[p+5] = 1;
}
inline void func61( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+2] = 1;
   seed[p+3] = 1;
   seed[p+4] = 1;
   seed[p+5] = 1;
   seed[p+7] = 1;
}
inline void func62( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+2] = 1;
   seed[p+3] = 1;
   seed[p+4] = 1;
   seed[p+5] = 1;
   seed[p+6] = 1;
}
inline void func63( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+2] = 1;
   seed[p+3] = 1;
   seed[p+4] = 1;
   seed[p+5] = 1;
   seed[p+6] = 1;
   seed[p+7] = 1;
}
inline void func64( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+1] = 1;
}
inline void func65( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+1] = 1;
   seed[p+7] = 1;
}
inline void func66( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+1] = 1;
   seed[p+6] = 1;
}
inline void func67( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+1] = 1;
   seed[p+6] = 1;
   seed[p+7] = 1;
}
inline void func68( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+1] = 1;
   seed[p+5] = 1;
}
inline void func69( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+1] = 1;
   seed[p+5] = 1;
   seed[p+7] = 1;
}
inline void func70( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+1] = 1;
   seed[p+5] = 1;
   seed[p+6] = 1;
}
inline void func71( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+1] = 1;
   seed[p+5] = 1;
   seed[p+6] = 1;
   seed[p+7] = 1;
}
inline void func72( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+1] = 1;
   seed[p+4] = 1;
}
inline void func73( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+1] = 1;
   seed[p+4] = 1;
   seed[p+7] = 1;
}
inline void func74( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+1] = 1;
   seed[p+4] = 1;
   seed[p+6] = 1;
}
inline void func75( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+1] = 1;
   seed[p+4] = 1;
   seed[p+6] = 1;
   seed[p+7] = 1;
}
inline void func76( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+1] = 1;
   seed[p+4] = 1;
   seed[p+5] = 1;
}
inline void func77( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+1] = 1;
   seed[p+4] = 1;
   seed[p+5] = 1;
   seed[p+7] = 1;
}
inline void func78( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+1] = 1;
   seed[p+4] = 1;
   seed[p+5] = 1;
   seed[p+6] = 1;
}
inline void func79( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+1] = 1;
   seed[p+4] = 1;
   seed[p+5] = 1;
   seed[p+6] = 1;
   seed[p+7] = 1;
}
inline void func80( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+1] = 1;
   seed[p+3] = 1;
}
inline void func81( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+1] = 1;
   seed[p+3] = 1;
   seed[p+7] = 1;
}
inline void func82( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+1] = 1;
   seed[p+3] = 1;
   seed[p+6] = 1;
}
inline void func83( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+1] = 1;
   seed[p+3] = 1;
   seed[p+6] = 1;
   seed[p+7] = 1;
}
inline void func84( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+1] = 1;
   seed[p+3] = 1;
   seed[p+5] = 1;
}
inline void func85( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+1] = 1;
   seed[p+3] = 1;
   seed[p+5] = 1;
   seed[p+7] = 1;
}
inline void func86( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+1] = 1;
   seed[p+3] = 1;
   seed[p+5] = 1;
   seed[p+6] = 1;
}
inline void func87( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+1] = 1;
   seed[p+3] = 1;
   seed[p+5] = 1;
   seed[p+6] = 1;
   seed[p+7] = 1;
}
inline void func88( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+1] = 1;
   seed[p+3] = 1;
   seed[p+4] = 1;
}
inline void func89( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+1] = 1;
   seed[p+3] = 1;
   seed[p+4] = 1;
   seed[p+7] = 1;
}
inline void func90( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+1] = 1;
   seed[p+3] = 1;
   seed[p+4] = 1;
   seed[p+6] = 1;
}
inline void func91( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+1] = 1;
   seed[p+3] = 1;
   seed[p+4] = 1;
   seed[p+6] = 1;
   seed[p+7] = 1;
}
inline void func92( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+1] = 1;
   seed[p+3] = 1;
   seed[p+4] = 1;
   seed[p+5] = 1;
}
inline void func93( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+1] = 1;
   seed[p+3] = 1;
   seed[p+4] = 1;
   seed[p+5] = 1;
   seed[p+7] = 1;
}
inline void func94( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+1] = 1;
   seed[p+3] = 1;
   seed[p+4] = 1;
   seed[p+5] = 1;
   seed[p+6] = 1;
}
inline void func95( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+1] = 1;
   seed[p+3] = 1;
   seed[p+4] = 1;
   seed[p+5] = 1;
   seed[p+6] = 1;
   seed[p+7] = 1;
}
inline void func96( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+1] = 1;
   seed[p+2] = 1;
}
inline void func97( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+1] = 1;
   seed[p+2] = 1;
   seed[p+7] = 1;
}
inline void func98( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+1] = 1;
   seed[p+2] = 1;
   seed[p+6] = 1;
}
inline void func99( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+1] = 1;
   seed[p+2] = 1;
   seed[p+6] = 1;
   seed[p+7] = 1;
}
inline void func100( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+1] = 1;
   seed[p+2] = 1;
   seed[p+5] = 1;
}
inline void func101( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+1] = 1;
   seed[p+2] = 1;
   seed[p+5] = 1;
   seed[p+7] = 1;
}
inline void func102( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+1] = 1;
   seed[p+2] = 1;
   seed[p+5] = 1;
   seed[p+6] = 1;
}
inline void func103( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+1] = 1;
   seed[p+2] = 1;
   seed[p+5] = 1;
   seed[p+6] = 1;
   seed[p+7] = 1;
}
inline void func104( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+1] = 1;
   seed[p+2] = 1;
   seed[p+4] = 1;
}
inline void func105( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+1] = 1;
   seed[p+2] = 1;
   seed[p+4] = 1;
   seed[p+7] = 1;
}
inline void func106( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+1] = 1;
   seed[p+2] = 1;
   seed[p+4] = 1;
   seed[p+6] = 1;
}
inline void func107( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+1] = 1;
   seed[p+2] = 1;
   seed[p+4] = 1;
   seed[p+6] = 1;
   seed[p+7] = 1;
}
inline void func108( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+1] = 1;
   seed[p+2] = 1;
   seed[p+4] = 1;
   seed[p+5] = 1;
}
inline void func109( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+1] = 1;
   seed[p+2] = 1;
   seed[p+4] = 1;
   seed[p+5] = 1;
   seed[p+7] = 1;
}
inline void func110( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+1] = 1;
   seed[p+2] = 1;
   seed[p+4] = 1;
   seed[p+5] = 1;
   seed[p+6] = 1;
}
inline void func111( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+1] = 1;
   seed[p+2] = 1;
   seed[p+4] = 1;
   seed[p+5] = 1;
   seed[p+6] = 1;
   seed[p+7] = 1;
}
inline void func112( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+1] = 1;
   seed[p+2] = 1;
   seed[p+3] = 1;
}
inline void func113( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+1] = 1;
   seed[p+2] = 1;
   seed[p+3] = 1;
   seed[p+7] = 1;
}
inline void func114( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+1] = 1;
   seed[p+2] = 1;
   seed[p+3] = 1;
   seed[p+6] = 1;
}
inline void func115( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+1] = 1;
   seed[p+2] = 1;
   seed[p+3] = 1;
   seed[p+6] = 1;
   seed[p+7] = 1;
}
inline void func116( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+1] = 1;
   seed[p+2] = 1;
   seed[p+3] = 1;
   seed[p+5] = 1;
}
inline void func117( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+1] = 1;
   seed[p+2] = 1;
   seed[p+3] = 1;
   seed[p+5] = 1;
   seed[p+7] = 1;
}
inline void func118( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+1] = 1;
   seed[p+2] = 1;
   seed[p+3] = 1;
   seed[p+5] = 1;
   seed[p+6] = 1;
}
inline void func119( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+1] = 1;
   seed[p+2] = 1;
   seed[p+3] = 1;
   seed[p+5] = 1;
   seed[p+6] = 1;
   seed[p+7] = 1;
}
inline void func120( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+1] = 1;
   seed[p+2] = 1;
   seed[p+3] = 1;
   seed[p+4] = 1;
}
inline void func121( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+1] = 1;
   seed[p+2] = 1;
   seed[p+3] = 1;
   seed[p+4] = 1;
   seed[p+7] = 1;
}
inline void func122( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+1] = 1;
   seed[p+2] = 1;
   seed[p+3] = 1;
   seed[p+4] = 1;
   seed[p+6] = 1;
}
inline void func123( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+1] = 1;
   seed[p+2] = 1;
   seed[p+3] = 1;
   seed[p+4] = 1;
   seed[p+6] = 1;
   seed[p+7] = 1;
}
inline void func124( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+1] = 1;
   seed[p+2] = 1;
   seed[p+3] = 1;
   seed[p+4] = 1;
   seed[p+5] = 1;
}
inline void func125( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+1] = 1;
   seed[p+2] = 1;
   seed[p+3] = 1;
   seed[p+4] = 1;
   seed[p+5] = 1;
   seed[p+7] = 1;
}
inline void func126( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+1] = 1;
   seed[p+2] = 1;
   seed[p+3] = 1;
   seed[p+4] = 1;
   seed[p+5] = 1;
   seed[p+6] = 1;
}
inline void func127( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+1] = 1;
   seed[p+2] = 1;
   seed[p+3] = 1;
   seed[p+4] = 1;
   seed[p+5] = 1;
   seed[p+6] = 1;
   seed[p+7] = 1;
}
inline void func128( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
}
inline void func129( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+7] = 1;
}
inline void func130( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+6] = 1;
}
inline void func131( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+6] = 1;
   seed[p+7] = 1;
}
inline void func132( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+5] = 1;
}
inline void func133( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+5] = 1;
   seed[p+7] = 1;
}
inline void func134( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+5] = 1;
   seed[p+6] = 1;
}
inline void func135( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+5] = 1;
   seed[p+6] = 1;
   seed[p+7] = 1;
}
inline void func136( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+4] = 1;
}
inline void func137( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+4] = 1;
   seed[p+7] = 1;
}
inline void func138( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+4] = 1;
   seed[p+6] = 1;
}
inline void func139( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+4] = 1;
   seed[p+6] = 1;
   seed[p+7] = 1;
}
inline void func140( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+4] = 1;
   seed[p+5] = 1;
}
inline void func141( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+4] = 1;
   seed[p+5] = 1;
   seed[p+7] = 1;
}
inline void func142( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+4] = 1;
   seed[p+5] = 1;
   seed[p+6] = 1;
}
inline void func143( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+4] = 1;
   seed[p+5] = 1;
   seed[p+6] = 1;
   seed[p+7] = 1;
}
inline void func144( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+3] = 1;
}
inline void func145( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+3] = 1;
   seed[p+7] = 1;
}
inline void func146( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+3] = 1;
   seed[p+6] = 1;
}
inline void func147( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+3] = 1;
   seed[p+6] = 1;
   seed[p+7] = 1;
}
inline void func148( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+3] = 1;
   seed[p+5] = 1;
}
inline void func149( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+3] = 1;
   seed[p+5] = 1;
   seed[p+7] = 1;
}
inline void func150( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+3] = 1;
   seed[p+5] = 1;
   seed[p+6] = 1;
}
inline void func151( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+3] = 1;
   seed[p+5] = 1;
   seed[p+6] = 1;
   seed[p+7] = 1;
}
inline void func152( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+3] = 1;
   seed[p+4] = 1;
}
inline void func153( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+3] = 1;
   seed[p+4] = 1;
   seed[p+7] = 1;
}
inline void func154( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+3] = 1;
   seed[p+4] = 1;
   seed[p+6] = 1;
}
inline void func155( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+3] = 1;
   seed[p+4] = 1;
   seed[p+6] = 1;
   seed[p+7] = 1;
}
inline void func156( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+3] = 1;
   seed[p+4] = 1;
   seed[p+5] = 1;
}
inline void func157( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+3] = 1;
   seed[p+4] = 1;
   seed[p+5] = 1;
   seed[p+7] = 1;
}
inline void func158( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+3] = 1;
   seed[p+4] = 1;
   seed[p+5] = 1;
   seed[p+6] = 1;
}
inline void func159( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+3] = 1;
   seed[p+4] = 1;
   seed[p+5] = 1;
   seed[p+6] = 1;
   seed[p+7] = 1;
}
inline void func160( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+2] = 1;
}
inline void func161( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+2] = 1;
   seed[p+7] = 1;
}
inline void func162( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+2] = 1;
   seed[p+6] = 1;
}
inline void func163( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+2] = 1;
   seed[p+6] = 1;
   seed[p+7] = 1;
}
inline void func164( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+2] = 1;
   seed[p+5] = 1;
}
inline void func165( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+2] = 1;
   seed[p+5] = 1;
   seed[p+7] = 1;
}
inline void func166( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+2] = 1;
   seed[p+5] = 1;
   seed[p+6] = 1;
}
inline void func167( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+2] = 1;
   seed[p+5] = 1;
   seed[p+6] = 1;
   seed[p+7] = 1;
}
inline void func168( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+2] = 1;
   seed[p+4] = 1;
}
inline void func169( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+2] = 1;
   seed[p+4] = 1;
   seed[p+7] = 1;
}
inline void func170( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+2] = 1;
   seed[p+4] = 1;
   seed[p+6] = 1;
}
inline void func171( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+2] = 1;
   seed[p+4] = 1;
   seed[p+6] = 1;
   seed[p+7] = 1;
}
inline void func172( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+2] = 1;
   seed[p+4] = 1;
   seed[p+5] = 1;
}
inline void func173( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+2] = 1;
   seed[p+4] = 1;
   seed[p+5] = 1;
   seed[p+7] = 1;
}
inline void func174( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+2] = 1;
   seed[p+4] = 1;
   seed[p+5] = 1;
   seed[p+6] = 1;
}
inline void func175( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+2] = 1;
   seed[p+4] = 1;
   seed[p+5] = 1;
   seed[p+6] = 1;
   seed[p+7] = 1;
}
inline void func176( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+2] = 1;
   seed[p+3] = 1;
}
inline void func177( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+2] = 1;
   seed[p+3] = 1;
   seed[p+7] = 1;
}
inline void func178( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+2] = 1;
   seed[p+3] = 1;
   seed[p+6] = 1;
}
inline void func179( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+2] = 1;
   seed[p+3] = 1;
   seed[p+6] = 1;
   seed[p+7] = 1;
}
inline void func180( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+2] = 1;
   seed[p+3] = 1;
   seed[p+5] = 1;
}
inline void func181( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+2] = 1;
   seed[p+3] = 1;
   seed[p+5] = 1;
   seed[p+7] = 1;
}
inline void func182( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+2] = 1;
   seed[p+3] = 1;
   seed[p+5] = 1;
   seed[p+6] = 1;
}
inline void func183( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+2] = 1;
   seed[p+3] = 1;
   seed[p+5] = 1;
   seed[p+6] = 1;
   seed[p+7] = 1;
}
inline void func184( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+2] = 1;
   seed[p+3] = 1;
   seed[p+4] = 1;
}
inline void func185( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+2] = 1;
   seed[p+3] = 1;
   seed[p+4] = 1;
   seed[p+7] = 1;
}
inline void func186( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+2] = 1;
   seed[p+3] = 1;
   seed[p+4] = 1;
   seed[p+6] = 1;
}
inline void func187( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+2] = 1;
   seed[p+3] = 1;
   seed[p+4] = 1;
   seed[p+6] = 1;
   seed[p+7] = 1;
}
inline void func188( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+2] = 1;
   seed[p+3] = 1;
   seed[p+4] = 1;
   seed[p+5] = 1;
}
inline void func189( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+2] = 1;
   seed[p+3] = 1;
   seed[p+4] = 1;
   seed[p+5] = 1;
   seed[p+7] = 1;
}
inline void func190( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+2] = 1;
   seed[p+3] = 1;
   seed[p+4] = 1;
   seed[p+5] = 1;
   seed[p+6] = 1;
}
inline void func191( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+2] = 1;
   seed[p+3] = 1;
   seed[p+4] = 1;
   seed[p+5] = 1;
   seed[p+6] = 1;
   seed[p+7] = 1;
}
inline void func192( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+1] = 1;
}
inline void func193( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+1] = 1;
   seed[p+7] = 1;
}
inline void func194( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+1] = 1;
   seed[p+6] = 1;
}
inline void func195( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+1] = 1;
   seed[p+6] = 1;
   seed[p+7] = 1;
}
inline void func196( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+1] = 1;
   seed[p+5] = 1;
}
inline void func197( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+1] = 1;
   seed[p+5] = 1;
   seed[p+7] = 1;
}
inline void func198( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+1] = 1;
   seed[p+5] = 1;
   seed[p+6] = 1;
}
inline void func199( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+1] = 1;
   seed[p+5] = 1;
   seed[p+6] = 1;
   seed[p+7] = 1;
}
inline void func200( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+1] = 1;
   seed[p+4] = 1;
}
inline void func201( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+1] = 1;
   seed[p+4] = 1;
   seed[p+7] = 1;
}
inline void func202( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+1] = 1;
   seed[p+4] = 1;
   seed[p+6] = 1;
}
inline void func203( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+1] = 1;
   seed[p+4] = 1;
   seed[p+6] = 1;
   seed[p+7] = 1;
}
inline void func204( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+1] = 1;
   seed[p+4] = 1;
   seed[p+5] = 1;
}
inline void func205( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+1] = 1;
   seed[p+4] = 1;
   seed[p+5] = 1;
   seed[p+7] = 1;
}
inline void func206( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+1] = 1;
   seed[p+4] = 1;
   seed[p+5] = 1;
   seed[p+6] = 1;
}
inline void func207( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+1] = 1;
   seed[p+4] = 1;
   seed[p+5] = 1;
   seed[p+6] = 1;
   seed[p+7] = 1;
}
inline void func208( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+1] = 1;
   seed[p+3] = 1;
}
inline void func209( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+1] = 1;
   seed[p+3] = 1;
   seed[p+7] = 1;
}
inline void func210( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+1] = 1;
   seed[p+3] = 1;
   seed[p+6] = 1;
}
inline void func211( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+1] = 1;
   seed[p+3] = 1;
   seed[p+6] = 1;
   seed[p+7] = 1;
}
inline void func212( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+1] = 1;
   seed[p+3] = 1;
   seed[p+5] = 1;
}
inline void func213( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+1] = 1;
   seed[p+3] = 1;
   seed[p+5] = 1;
   seed[p+7] = 1;
}
inline void func214( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+1] = 1;
   seed[p+3] = 1;
   seed[p+5] = 1;
   seed[p+6] = 1;
}
inline void func215( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+1] = 1;
   seed[p+3] = 1;
   seed[p+5] = 1;
   seed[p+6] = 1;
   seed[p+7] = 1;
}
inline void func216( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+1] = 1;
   seed[p+3] = 1;
   seed[p+4] = 1;
}
inline void func217( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+1] = 1;
   seed[p+3] = 1;
   seed[p+4] = 1;
   seed[p+7] = 1;
}
inline void func218( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+1] = 1;
   seed[p+3] = 1;
   seed[p+4] = 1;
   seed[p+6] = 1;
}
inline void func219( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+1] = 1;
   seed[p+3] = 1;
   seed[p+4] = 1;
   seed[p+6] = 1;
   seed[p+7] = 1;
}
inline void func220( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+1] = 1;
   seed[p+3] = 1;
   seed[p+4] = 1;
   seed[p+5] = 1;
}
inline void func221( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+1] = 1;
   seed[p+3] = 1;
   seed[p+4] = 1;
   seed[p+5] = 1;
   seed[p+7] = 1;
}
inline void func222( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+1] = 1;
   seed[p+3] = 1;
   seed[p+4] = 1;
   seed[p+5] = 1;
   seed[p+6] = 1;
}
inline void func223( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+1] = 1;
   seed[p+3] = 1;
   seed[p+4] = 1;
   seed[p+5] = 1;
   seed[p+6] = 1;
   seed[p+7] = 1;
}
inline void func224( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+1] = 1;
   seed[p+2] = 1;
}
inline void func225( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+1] = 1;
   seed[p+2] = 1;
   seed[p+7] = 1;
}
inline void func226( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+1] = 1;
   seed[p+2] = 1;
   seed[p+6] = 1;
}
inline void func227( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+1] = 1;
   seed[p+2] = 1;
   seed[p+6] = 1;
   seed[p+7] = 1;
}
inline void func228( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+1] = 1;
   seed[p+2] = 1;
   seed[p+5] = 1;
}
inline void func229( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+1] = 1;
   seed[p+2] = 1;
   seed[p+5] = 1;
   seed[p+7] = 1;
}
inline void func230( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+1] = 1;
   seed[p+2] = 1;
   seed[p+5] = 1;
   seed[p+6] = 1;
}
inline void func231( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+1] = 1;
   seed[p+2] = 1;
   seed[p+5] = 1;
   seed[p+6] = 1;
   seed[p+7] = 1;
}
inline void func232( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+1] = 1;
   seed[p+2] = 1;
   seed[p+4] = 1;
}
inline void func233( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+1] = 1;
   seed[p+2] = 1;
   seed[p+4] = 1;
   seed[p+7] = 1;
}
inline void func234( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+1] = 1;
   seed[p+2] = 1;
   seed[p+4] = 1;
   seed[p+6] = 1;
}
inline void func235( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+1] = 1;
   seed[p+2] = 1;
   seed[p+4] = 1;
   seed[p+6] = 1;
   seed[p+7] = 1;
}
inline void func236( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+1] = 1;
   seed[p+2] = 1;
   seed[p+4] = 1;
   seed[p+5] = 1;
}
inline void func237( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+1] = 1;
   seed[p+2] = 1;
   seed[p+4] = 1;
   seed[p+5] = 1;
   seed[p+7] = 1;
}
inline void func238( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+1] = 1;
   seed[p+2] = 1;
   seed[p+4] = 1;
   seed[p+5] = 1;
   seed[p+6] = 1;
}
inline void func239( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+1] = 1;
   seed[p+2] = 1;
   seed[p+4] = 1;
   seed[p+5] = 1;
   seed[p+6] = 1;
   seed[p+7] = 1;
}
inline void func240( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+1] = 1;
   seed[p+2] = 1;
   seed[p+3] = 1;
}
inline void func241( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+1] = 1;
   seed[p+2] = 1;
   seed[p+3] = 1;
   seed[p+7] = 1;
}
inline void func242( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+1] = 1;
   seed[p+2] = 1;
   seed[p+3] = 1;
   seed[p+6] = 1;
}
inline void func243( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+1] = 1;
   seed[p+2] = 1;
   seed[p+3] = 1;
   seed[p+6] = 1;
   seed[p+7] = 1;
}
inline void func244( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+1] = 1;
   seed[p+2] = 1;
   seed[p+3] = 1;
   seed[p+5] = 1;
}
inline void func245( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+1] = 1;
   seed[p+2] = 1;
   seed[p+3] = 1;
   seed[p+5] = 1;
   seed[p+7] = 1;
}
inline void func246( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+1] = 1;
   seed[p+2] = 1;
   seed[p+3] = 1;
   seed[p+5] = 1;
   seed[p+6] = 1;
}
inline void func247( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+1] = 1;
   seed[p+2] = 1;
   seed[p+3] = 1;
   seed[p+5] = 1;
   seed[p+6] = 1;
   seed[p+7] = 1;
}
inline void func248( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+1] = 1;
   seed[p+2] = 1;
   seed[p+3] = 1;
   seed[p+4] = 1;
}
inline void func249( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+1] = 1;
   seed[p+2] = 1;
   seed[p+3] = 1;
   seed[p+4] = 1;
   seed[p+7] = 1;
}
inline void func250( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+1] = 1;
   seed[p+2] = 1;
   seed[p+3] = 1;
   seed[p+4] = 1;
   seed[p+6] = 1;
}
inline void func251( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+1] = 1;
   seed[p+2] = 1;
   seed[p+3] = 1;
   seed[p+4] = 1;
   seed[p+6] = 1;
   seed[p+7] = 1;
}
inline void func252( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+1] = 1;
   seed[p+2] = 1;
   seed[p+3] = 1;
   seed[p+4] = 1;
   seed[p+5] = 1;
}
inline void func253( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+1] = 1;
   seed[p+2] = 1;
   seed[p+3] = 1;
   seed[p+4] = 1;
   seed[p+5] = 1;
   seed[p+7] = 1;
}
inline void func254( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+1] = 1;
   seed[p+2] = 1;
   seed[p+3] = 1;
   seed[p+4] = 1;
   seed[p+5] = 1;
   seed[p+6] = 1;
}
inline void func255( int* seed, int i, int j ){
   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8
   seed[p+0] = 1;
   seed[p+1] = 1;
   seed[p+2] = 1;
   seed[p+3] = 1;
   seed[p+4] = 1;
   seed[p+5] = 1;
   seed[p+6] = 1;
   seed[p+7] = 1;
}



void (*find_seed_pos [])( int *seed,int i,int j) = {
   func0,func1,func2,func3,func4,func5,func6,func7,func8,func9,func10,func11,func12,func13,func14,func15,
   func16,func17,func18,func19,func20,func21,func22,func23,func24,func25,func26,func27,func28,func29,func30,func31,
   func32,func33,func34,func35,func36,func37,func38,func39,func40,func41,func42,func43,func44,func45,func46,func47,
   func48,func49,func50,func51,func52,func53,func54,func55,func56,func57,func58,func59,func60,func61,func62,func63,
   func64,func65,func66,func67,func68,func69,func70,func71,func72,func73,func74,func75,func76,func77,func78,func79,
   func80,func81,func82,func83,func84,func85,func86,func87,func88,func89,func90,func91,func92,func93,func94,func95,
   func96,func97,func98,func99,func100,func101,func102,func103,func104,func105,func106,func107,func108,func109,func110,func111,
   func112,func113,func114,func115,func116,func117,func118,func119,func120,func121,func122,func123,func124,func125,func126,func127,
   func128,func129,func130,func131,func132,func133,func134,func135,func136,func137,func138,func139,func140,func141,func142,func143,
   func144,func145,func146,func147,func148,func149,func150,func151,func152,func153,func154,func155,func156,func157,func158,func159,
   func160,func161,func162,func163,func164,func165,func166,func167,func168,func169,func170,func171,func172,func173,func174,func175,
   func176,func177,func178,func179,func180,func181,func182,func183,func184,func185,func186,func187,func188,func189,func190,func191,
   func192,func193,func194,func195,func196,func197,func198,func199,func200,func201,func202,func203,func204,func205,func206,func207,
   func208,func209,func210,func211,func212,func213,func214,func215,func216,func217,func218,func219,func220,func221,func222,func223,
   func224,func225,func226,func227,func228,func229,func230,func231,func232,func233,func234,func235,func236,func237,func238,func239,
   func240,func241,func242,func243,func244,func245,func246,func247,func248,func249,func250,func251,func252,func253,func254,func255
};
