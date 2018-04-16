// Include files for SSE4
#include "emmintrin.h"
#include "xmmintrin.h"
#include "immintrin.h"
#include "tmmintrin.h"
#include "emmintrin.h"
#include "time_meas.h"

#include <stdint.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>

#define VECTOR_SIZE 128 // TODO check this
#define SPARE 32 // TODO check this

// example SIMD macros, not necessary to be used, write your own

extern __m128i *reflip;
inline void cmult(__m128i a,__m128i b, __m128i *re32, __m128i *im32) {

    __m128i mmtmpb;

    mmtmpb    = _mm_sign_epi16(b,*(__m128i*)reflip); 
    *re32     = _mm_madd_epi16(a,mmtmpb);
    mmtmpb    = _mm_shufflelo_epi16(b,_MM_SHUFFLE(2,3,0,1));
    mmtmpb    = _mm_shufflehi_epi16(mmtmpb,_MM_SHUFFLE(2,3,0,1));
    *im32  = _mm_madd_epi16(a,mmtmpb);

}

inline __m128i cpack(__m128i xre,__m128i xim) {

    __m128i cpack_tmp1, cpack_tmp2;
    cpack_tmp1 = _mm_unpacklo_epi32(xre,xim);
    cpack_tmp1 = _mm_srai_epi32(cpack_tmp1,15);
    cpack_tmp2 = _mm_unpackhi_epi32(xre,xim);
    cpack_tmp2 = _mm_srai_epi32(cpack_tmp2,15);
    return(_mm_packs_epi32(cpack_tmp1,cpack_tmp2));

}

// routines to be written
void componentwise_multiply_real_scalar(int16_t *x,int16_t *y,int16_t *z,uint16_t N) {
    int16_t scal_tmp[8];

    for (int i = 0; i < N; i++) {
        z[i] = x[i] * y[i];
    }
}

// routines to be written
void componentwise_multiply_real_sse4(int16_t *x,int16_t *y,int16_t *z,uint16_t N) {

    __m128i *x128 = (__m128i *)x;
    __m128i *y128 = (__m128i *)y;
    __m128i *z128 = (__m128i *)z;

    for (int i = 0; i < N/8; i++) {
       z128[i] = _mm_mulhrs_epi16(x128[i], y128[i]); 
    }

}

int main() {
    int16_t *x, *y, *z;
    /*__m128i *x128, *y128, *z128;*/

    time_stats_t *t;

    x = malloc( (VECTOR_SIZE+SPARE) * sizeof(int16_t));
    y = malloc( (VECTOR_SIZE+SPARE) * sizeof(int16_t));
    z = malloc( (VECTOR_SIZE+SPARE) * sizeof(int16_t));

    // FILL THE ARRAY WITH RANDOM NUMBERS
    srand(time(NULL)); 
    for (int i = 0; i < VECTOR_SIZE + SPARE; i++) {
        x[i] = (int16_t)rand();
        y[i] = (int16_t)rand();
        /*printf("x: %d y: %d\n", x[i], y[i]);*/
    }
    

    /*x128 = malloc( (int)((VECTOR_SIZE+SPARE)/4) * sizeof(__m128i));*/
    /*y128 = malloc( (int)((VECTOR_SIZE+SPARE)/4) * sizeof(__m128i));*/
    /*z128 = malloc( (int)((VECTOR_SIZE+SPARE)/4) * sizeof(__m128i));*/

    // FOR testing
    /*x[0] = 1<<15;*/
    /*y[0] = 1;*/


    start_meas(t);
    componentwise_multiply_real_scalar(x, y, z, VECTOR_SIZE);
    stop_meas(t);

    reset_meas(t);

    // FOR DEBUGGING
    /*for (int i = 0; i < VECTOR_SIZE+SPARE; i++) {*/
        /*printf("It %d: %d * %d = %d\n", i, x[i], y[i], z[i]);*/
    /*}*/

    start_meas(t);
    componentwise_multiply_real_sse4(x, y, z, VECTOR_SIZE);
    stop_meas(t);
    
    return 0;
}
