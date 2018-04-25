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
#include <limits.h>
#include <math.h>
#include <stdlib.h>

#define VECTOR_SIZE 50000 // TODO check this
#define SPARE 288 // TODO check this
#define ITER_NUM 50

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
void componentwise_multiply_real_scalar(int16_t *x,int16_t *y,int16_t *z, int N) {

	int i;
	for (i = 0; i < N; i++) {
		z[i] = x[i] * y[i];
	}
}

// 128 bits inline function
void componentwise_multiply_real_sse4(int16_t *x,int16_t *y,int16_t *z, int N) {

	__m128i *x128 = (__m128i *)x;
	__m128i *y128 = (__m128i *)y;
	__m128i *z128 = (__m128i *)z;

	int i;
	for (i = 0; i < ceil(N/8.0); i+=1) {
		z128[i] = _mm_mulhrs_epi16(x128[i], y128[i]); 
	}

}

// 256 bits inline function
void componentwise_multiply_real_avx2(int16_t *x,int16_t *y,int16_t *z, uint16_t N) {

	__m256i *x256 = (__m256i *)x;
	__m256i *y256 = (__m256i *)y;
	__m256i *z256 = (__m256i *)z;

	int i;
	for (i = 0; i < ceil(N/16.0); i+=1) {
		z256[i] = _mm256_mulhrs_epi16(x256[i], y256[i]);
	}

}

int main() {
	/*int16_t *x, *y, *z;*/
	int16_t *x __attribute__((aligned(32)));
	int16_t *y __attribute__((aligned(32)));
	int16_t *z __attribute__((aligned(32)));
	FILE *f;
	/*__m128i *x128, *y128, *z128;*/

	time_stats_t t;

	f = fopen("measurements", "w");
	x = (int16_t *)malloc( (VECTOR_SIZE+SPARE) * sizeof(int16_t));
	y = (int16_t *)malloc( (VECTOR_SIZE+SPARE) * sizeof(int16_t));
	z = (int16_t *)malloc( (VECTOR_SIZE+SPARE) * sizeof(int16_t));

	/*x = (int16_t *) aligned_alloc(32,sizeof(int16_t) * VECTOR_SIZE);*/
	/*y = (int16_t *) aligned_alloc(32,sizeof(int16_t) * VECTOR_SIZE);*/
	/*z = (int16_t *) aligned_alloc(32,sizeof(int16_t) * VECTOR_SIZE);*/

	// FILL THE ARRAY WITH RANDOM NUMBERS
	int i;
	srand(time(NULL)); 
	for (i = 0; i < VECTOR_SIZE + SPARE; i++) {
		x[i] = (int16_t)rand();
		y[i] = (int16_t)rand();
		/*printf("x: %d y: %d\n", x[i], y[i]);*/
	}


	// FOR testing
	/*x[0] = 1<<15;*/
	/*y[0] = 1;*/


	long long int min;
	int j = 0;
	for (i = 0; i < VECTOR_SIZE; ++i) {
		min = INT_MAX;
		for (j = 0; j < ITER_NUM; ++j) {
			start_meas(&t);
			componentwise_multiply_real_scalar(x, y, z, i);
			stop_meas(&t);
			if (min > t.diff) {
				min = t.diff;
			}
			reset_meas(&t);
		}
		fprintf(f, "scalar %lld %d\n", min, i); 
	}
	printf("DONE SCALAR\n");


	// FOR DEBUGGING
	/*for (int i = 0; i < VECTOR_SIZE+SPARE; i++) {*/
	/*printf("It %d: %d * %d = %d\n", i, x[i], y[i], z[i]);*/
	/*}*/
	
	// 128 bits Multiplication
	for (i = 0; i < VECTOR_SIZE; ++i) {
		min = LLONG_MAX;
		for (j = 0; j < ITER_NUM; ++j) {
			start_meas(&t);
			componentwise_multiply_real_sse4(x, y, z, i);
			stop_meas(&t);
			if (min > t.diff) {
				min = t.diff;
			}
			reset_meas(&t);
		}
		fprintf(f, "vector %lld %d\n", min, i); 
	}
	printf("DONE 128\n");

	// 256 bits Multiplication
	for (i = 0; i < VECTOR_SIZE; ++i) {
		min = LLONG_MAX;
		for (j = 0; j < ITER_NUM; ++j) {
			start_meas(&t);
			componentwise_multiply_real_avx2(x, y, z, i);
			stop_meas(&t);
			if (min > t.diff) {
				min = t.diff;
			}
			reset_meas(&t);
		}
		fprintf(f, "vector256 %lld %d\n", min, i); 
	}
	printf("DONE 256\n");
	return 0;
}
