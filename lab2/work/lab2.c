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
#include <string.h>

#define VECTOR_SIZE 10000 // TODO check this
#define ITER_NUM 1

// example SIMD macros, not necessary to be used, write your own
//

__m128i reflip[16]={{1, -1}, {1, -1}, {1, -1}, {1, -1}, {1, -1}};

typedef struct cstruct{
	int16_t re;
	int16_t im;
}cstruct;

static inline void cmult(__m128i a,__m128i b, __m128i *re32, __m128i *im32);
static inline __m128i cpack(__m128i xre,__m128i xim);

static inline void cmult256(__m256i a,__m256i b, __m256i *re32, __m256i *im32);
static inline __m256i cpack256(__m256i xre,__m256i xim);

static inline void componentwise_multiply_complex_scalar(cstruct *x, cstruct *y, cstruct *z, int N);
static inline void componentwise_multiply_complex_128(cstruct *x, cstruct *y, cstruct *z, int N);
static inline void componentwise_multiply_complex_256(cstruct *x, cstruct *y, cstruct *z, int N);

static inline void componentwise_multiply_complex_scalarB(int16_t *xre, int16_t *yre, int16_t *xim, int16_t *yim, int16_t *zre, int16_t *zim, int N);
static inline void componentwise_multiply_complex_128B(int16_t *xre, int16_t *yre, int16_t *xim, int16_t *yim, int16_t *zre, int16_t *zim, int N);
static inline void componentwise_multiply_complex_256B(int16_t *xre, int16_t *yre, int16_t *xim, int16_t *yim, int16_t *zre, int16_t *zim, int N);

static inline void componentwise_multiply_real_scalar(int16_t *x,int16_t *y,int16_t *z, int N);
static inline void componentwise_multiply_real_128(int16_t *x,int16_t *y,int16_t *z, int N);
static inline void componentwise_multiply_real_256(int16_t *x,int16_t *y,int16_t *z, uint16_t N);


int main(int argc, char *argv[]) {
	/*int16_t *x, *y, *z;*/
	time_stats_t t;
	long long int acc;
	long long int min;
	int i = 0;
	int j = 0;

	int16_t *x __attribute__((aligned(32)));
	int16_t *y __attribute__((aligned(32)));
	int16_t *z __attribute__((aligned(32)));
	int16_t *re1 __attribute__((aligned(32)));
	int16_t *re2 __attribute__((aligned(32)));
	int16_t *re3 __attribute__((aligned(32)));
	int16_t *im1 __attribute__((aligned(32)));
	int16_t *im2 __attribute__((aligned(32)));
	int16_t *im3 __attribute__((aligned(32)));
	cstruct *xc __attribute__((aligned(32)));
	cstruct *yc __attribute__((aligned(32))); 
	cstruct *zc __attribute__((aligned(32)));
	FILE *f;
	f = fopen(argv[2], "w");
	/*__m128i *x128, *y128, *z128;*/

	if (argc != 4) {
		fprintf(stderr, "Usage ./lab2 [avg/min] [filename] [real/complex/complexB]\n");
		return(-1);
	}

	if (!strcmp(argv[3], "complex")){
		printf("Complex selected\n");

		xc = (cstruct *) aligned_alloc(32,sizeof(cstruct) * VECTOR_SIZE);
		yc = (cstruct *) aligned_alloc(32,sizeof(cstruct) * VECTOR_SIZE);
		zc = (cstruct *) aligned_alloc(32,sizeof(cstruct) * VECTOR_SIZE);

		// FILL THE ARRAY WITH RANDOM NUMBERS
		int i;
		srand(time(NULL)); 
		for (i = 0; i < VECTOR_SIZE; i++) {
			xc[i].re = (int16_t)rand();
			xc[i].im = (int16_t)rand();
			yc[i].re = (int16_t)rand();
			yc[i].im = (int16_t)rand();
			/*printf("x: %d y: %d\n", x[i], y[i]);*/
		}
		printf("finished init\n");
	}
	else if (!strcmp(argv[3], "complexB")) {

		re1 = (int16_t *) aligned_alloc(32,sizeof(int16_t) * VECTOR_SIZE);
		re2 = (int16_t *) aligned_alloc(32,sizeof(int16_t) * VECTOR_SIZE);
		re3 = (int16_t *) aligned_alloc(32,sizeof(int16_t) * VECTOR_SIZE);
		im1 = (int16_t *) aligned_alloc(32,sizeof(int16_t) * VECTOR_SIZE);
		im2 = (int16_t *) aligned_alloc(32,sizeof(int16_t) * VECTOR_SIZE);
		im3 = (int16_t *) aligned_alloc(32,sizeof(int16_t) * VECTOR_SIZE);

		// FILL THE ARRAY WITH RANDOM NUMBERS
		int i;
		srand(time(NULL)); 
		for (i = 0; i < VECTOR_SIZE; i++) {
			re1[i] = (int16_t)rand();
			re2[i] = (int16_t)rand();
			im1[i] = (int16_t)rand();
			im2[i] = (int16_t)rand();
			/*printf("x: %d y: %d\n", x[i], y[i]);*/
		}

	}
	else {
		x = (int16_t *) aligned_alloc(32,sizeof(int16_t) * VECTOR_SIZE);
		y = (int16_t *) aligned_alloc(32,sizeof(int16_t) * VECTOR_SIZE);
		z = (int16_t *) aligned_alloc(32,sizeof(int16_t) * VECTOR_SIZE);
		// FILL THE ARRAY WITH RANDOM NUMBERS
		int i;
		srand(time(NULL)); 
		for (i = 0; i < VECTOR_SIZE; i++) {
			x[i] = (int16_t)rand();
			y[i] = (int16_t)rand();
			/*printf("x: %d y: %d\n", x[i], y[i]);*/
		}
	}


	// FOR testing
	/*x[0] = 1<<15;*/
	/*y[0] = 1;*/


	if (!strcmp(argv[3], "complex")) {
		for (i = 0; i < VECTOR_SIZE; ++i) {
			min = INT_MAX;
			acc = 0;
			for (j = 0; j < ITER_NUM; ++j) {
				start_meas(&t);
				componentwise_multiply_complex_scalar(xc, yc, zc, i);
				stop_meas(&t);
				if (min > t.diff) {
					min = t.diff;
				}
				acc+=t.diff;
				reset_meas(&t);
			}
			if (!strcmp(argv[1],"avg")) {
				printf("scalar %f %d\n", (float)(acc/ITER_NUM), i); 
			}
			else {
				printf("scalar %lld %d\n", min, i); 
			}
		}
		printf("DONE SCALAR\n");
		// 128 bits complex Multiplication
		for (i = 0; i < VECTOR_SIZE; ++i) {
			min = LLONG_MAX;
			acc = 0;
			for (j = 0; j < ITER_NUM; ++j) {
				start_meas(&t);
				componentwise_multiply_complex_128(xc, yc, zc, i);
				stop_meas(&t);
				if (min > t.diff) {
					min = t.diff;
				}
				acc+=t.diff;
				reset_meas(&t);
			}
			if (strcmp(argv[1],"avg")) {
				printf("vector %f %d\n", (float)(acc/ITER_NUM), i); 
			}
			else {
				printf("vector %lld %d\n", min, i); 
			}
		}
		printf("DONE 128\n");

		// 256 bits complex Multiplication
		for (i = 0; i < VECTOR_SIZE; ++i) {
			min = LLONG_MAX;
			acc = 0;
			for (j = 0; j < ITER_NUM; ++j) {
				start_meas(&t);
				componentwise_multiply_complex_256(xc, yc, zc, i);
				stop_meas(&t);
				if (min > t.diff) {
					min = t.diff;
				}
				acc+=t.diff;
				reset_meas(&t);
			}
			if (strcmp(argv[1],"avg")) {
				printf("vector256 %f %d\n", (float)(acc/ITER_NUM), i); 
			}
			else {
				printf("vector256 %lld %d\n", min, i); 
			}
		}
		printf("DONE 256\n");

	}
	else if (!strcmp(argv[3], "complexB")) {
		for (i = 0; i < VECTOR_SIZE; ++i) {
			min = INT_MAX;
			acc = 0;
			for (j = 0; j < ITER_NUM; ++j) {
				start_meas(&t);
				componentwise_multiply_complex_scalarB(re1, re2, im1, im2, re3, im3, i);
				stop_meas(&t);
				if (min > t.diff) {
					min = t.diff;
				}
				acc+=t.diff;
				reset_meas(&t);
			}
			if (!strcmp(argv[1],"avg")) {
				printf("scalar %f %d\n", (float)(acc/ITER_NUM), i); 
			}
			else {
				printf("scalar %lld %d\n", min, i); 
			}
		}
		printf("DONE SCALAR\n");
		// 128 bits complex Multiplication
		for (i = 0; i < VECTOR_SIZE; ++i) {
			min = LLONG_MAX;
			acc = 0;
			for (j = 0; j < ITER_NUM; ++j) {
				start_meas(&t);
				componentwise_multiply_complex_128B(re1, re2, im1, im2, re3, im3, i);
				stop_meas(&t);
				componentwise_multiply_complex_128B(re1, re2, im1, im2, re3, im3, i);
				if (min > t.diff) {
					min = t.diff;
				}
				acc+=t.diff;
				reset_meas(&t);
			}
			if (strcmp(argv[1],"avg")) {
				printf("vector %f %d\n", (float)(acc/ITER_NUM), i); 
			}
			else {
				printf("vector %lld %d\n", min, i); 
			}
		}
		printf("DONE 128\n");

		// 256 bits complex Multiplication
		for (i = 0; i < VECTOR_SIZE; ++i) {
			min = LLONG_MAX;
			acc = 0;
			for (j = 0; j < ITER_NUM; ++j) {
				start_meas(&t);
				componentwise_multiply_complex_256B(re1, re2, im1, im2, re3, im3, i);
				stop_meas(&t);
				if (min > t.diff) {
					min = t.diff;
				}
				acc+=t.diff;
				reset_meas(&t);
			}
			if (strcmp(argv[1],"avg")) {
				printf("vector256 %f %d\n", (float)(acc/ITER_NUM), i); 
			}
			else {
				printf("vector256 %lld %d\n", min, i); 
			}
		}
		printf("DONE 256\n");

	}
	else {
		for (i = 0; i < VECTOR_SIZE; ++i) {
			min = INT_MAX;
			acc = 0;
			for (j = 0; j < ITER_NUM; ++j) {
				start_meas(&t);
				componentwise_multiply_real_scalar(x, y, z, i);
				stop_meas(&t);
				if (min > t.diff) {
					min = t.diff;
				}
				acc+=t.diff;
				reset_meas(&t);
			}
			if (!strcmp(argv[1],"avg")) {
				printf("scalar %f %d\n", (float)(acc/ITER_NUM), i); 

			}
			else {
				printf("scalar %lld %d\n", min, i); 
			}
		}
		printf("DONE SCALAR\n");


		// FOR DEBUGGING
		/*for (int i = 0; i < VECTOR_SIZE+SPARE; i++) {*/
		/*printf("It %d: %d * %d = %d\n", i, x[i], y[i], z[i]);*/
		/*}*/

		// 128 bits Multiplication
		for (i = 0; i < VECTOR_SIZE; ++i) {
			min = LLONG_MAX;
			acc = 0;
			for (j = 0; j < ITER_NUM; ++j) {
				start_meas(&t);
				componentwise_multiply_real_128(x, y, z, i);
				stop_meas(&t);
				if (min > t.diff) {
					min = t.diff;
				}
				acc+=t.diff;
				reset_meas(&t);
			}
			if (strcmp(argv[1],"avg")) {
				printf("vector %f %d\n", (float)(acc/ITER_NUM), i); 
			}
			else {
				printf("vector %lld %d\n", min, i); 
			}
		}
		printf("DONE 128\n");

		// 256 bits Multiplication
		for (i = 0; i < VECTOR_SIZE; ++i) {
			min = LLONG_MAX;
			acc = 0;
			for (j = 0; j < ITER_NUM; ++j) {
				start_meas(&t);
				componentwise_multiply_real_256(x, y, z, i);
				stop_meas(&t);
				if (min > t.diff) {
					min = t.diff;
				}
				acc+=t.diff;
				reset_meas(&t);
			}
			if (strcmp(argv[1],"avg")) {
				printf("vector256 %f %d\n", (float)(acc/ITER_NUM), i); 
			}
			else {
				printf("vector256 %lld %d\n", min, i); 
			}
		}
		printf("DONE 256\n");
	}
	return 0;
}

static inline void cmult(__m128i a,__m128i b, __m128i *re32, __m128i *im32) {
	__m128i mmtmpb;
	/*printf("cmult\n");*/
	mmtmpb    = _mm_sign_epi16(b,*(__m128i*)reflip); 
	*re32     = _mm_madd_epi16(a,mmtmpb);
	mmtmpb    = _mm_shufflelo_epi16(b,_MM_SHUFFLE(2,3,0,1));
	mmtmpb    = _mm_shufflehi_epi16(mmtmpb,_MM_SHUFFLE(2,3,0,1));
	*im32  = _mm_madd_epi16(a,mmtmpb);
}

static inline __m128i cpack(__m128i xre, __m128i xim) {
	__m128i cpack_tmp1, cpack_tmp2;
	cpack_tmp1 = _mm_unpacklo_epi32(xre,xim);
	cpack_tmp1 = _mm_srai_epi32(cpack_tmp1,15);
	cpack_tmp2 = _mm_unpackhi_epi32(xre,xim);
	cpack_tmp2 = _mm_srai_epi32(cpack_tmp2,15);
	return(_mm_packs_epi32(cpack_tmp1,cpack_tmp2));
}

static inline void cmult256(__m256i a,__m256i b, __m256i *re32, __m256i *im32) {
	__m256i mmtmpb;
	/*printf("cmult\n");*/
	mmtmpb    = _mm256_sign_epi16(b,*(__m256i*)reflip); 
	*re32     = _mm256_madd_epi16(a,mmtmpb);
	mmtmpb    = _mm256_shufflelo_epi16(b,_MM_SHUFFLE(2,3,0,1));
	mmtmpb    = _mm256_shufflehi_epi16(mmtmpb,_MM_SHUFFLE(2,3,0,1));
	*im32  = _mm256_madd_epi16(a,mmtmpb);
}

static inline __m256i cpack256(__m256i xre, __m256i xim) {
	__m256i cpack_tmp1, cpack_tmp2;
	cpack_tmp1 = _mm256_unpacklo_epi32(xre,xim);
	cpack_tmp1 = _mm256_srai_epi32(cpack_tmp1,15);
	cpack_tmp2 = _mm256_unpackhi_epi32(xre,xim);
	cpack_tmp2 = _mm256_srai_epi32(cpack_tmp2,15);
	return(_mm256_packs_epi32(cpack_tmp1,cpack_tmp2));
}


static inline void componentwise_multiply_complex_scalarB(int16_t *xre, int16_t *yre, int16_t *xim, int16_t *yim, int16_t *zre, int16_t *zim, int N) {
	int i;
	for (i = 0; i < N; i++) {
		zre[i] = xre[i]*yre[i] - xim[i]*yim[i];
		zim[i] = xre[i]*yim[i] + xim[i]*yre[i];
	}
}

static inline void componentwise_multiply_complex_128B(int16_t *xre, int16_t *yre, int16_t *xim, int16_t *yim, int16_t *zre, int16_t *zim, int N) {
	__m128i *xre128 = (__m128i *)xre;
	__m128i *yre128 = (__m128i *)yre;
	__m128i *xim128 = (__m128i *)xim;
	__m128i *yim128 = (__m128i *)yim;
	__m128i *zre128 = (__m128i *)zre;
	__m128i *zim128 = (__m128i *)zim;
	int i;
	for (i = 0; i < N>>8; i++) {
		zre128[i] = _mm_mulhrs_epi16(xre128[i], yre128[i]); 
		zre128[i] = _mm_sub_epi16(zre128[i], _mm_mulhrs_epi16(xim128[i], yim128[i])); 
		zim128[i] = _mm_mulhrs_epi16(xre128[i], yim128[i]); 
		zim128[i] = _mm_add_epi16(zre128[i], _mm_mulhrs_epi16(xim128[i], yre128[i])); 
	}
}

static inline void componentwise_multiply_complex_256B(int16_t *xre, int16_t *yre, int16_t *xim, int16_t *yim, int16_t *zre, int16_t *zim, int N) {
	__m256i *xre256 = (__m256i *)xre;
	__m256i *yre256 = (__m256i *)yre;
	__m256i *xim256 = (__m256i *)xim;
	__m256i *yim256 = (__m256i *)yim;
	__m256i *zre256 = (__m256i *)zre;
	__m256i *zim256 = (__m256i *)zim;
	int i;
	for (i = 0; i < (N>>16); i++) {
		zre256[i] = _mm256_mulhrs_epi16(xre256[i], yre256[i]); 
		zre256[i] = _mm256_sub_epi16(zre256[i], _mm256_mulhrs_epi16(xim256[i], yim256[i])); 
		zim256[i] = _mm256_mulhrs_epi16(xre256[i], yim256[i]); 
		zim256[i] = _mm256_add_epi16(zre256[i], _mm256_mulhrs_epi16(xim256[i], yre256[i])); 
	}
}

static inline void componentwise_multiply_complex_scalar(cstruct *x, cstruct *y, cstruct *z, int N) {
	int i;
	for (i = 0; i < N; i++) {
		z[i].re = x[i].re * y[i].re - x[i].im * y[i].im;
		z[i].im = x[i].re * y[i].im + x[i].im * y[i].re;
	}
}

static inline void componentwise_multiply_complex_128(cstruct *x, cstruct *y, cstruct *z, int N) {
	__m128i *x128 = (__m128i *)x;
	__m128i *y128 = (__m128i *)y;
	__m128i *z128 = (__m128i *)z;
	__m128i res_re;
	__m128i res_im;
	int i;
	for (i = 0; i < ceil(N/8.0)*2; i++) {
		cmult(x128[i], y128[i], &res_re, &res_im); 
		z128[i] = cpack(res_re, res_im); 
	}
}

static inline void componentwise_multiply_complex_256(cstruct *x, cstruct *y, cstruct *z, int N) {
	__m256i *x256 = (__m256i *)x;
	__m256i *y256 = (__m256i *)y;
	__m256i *z256 = (__m256i *)z;
	__m256i res_re;
	__m256i res_im;
	int i;
	for (i = 0; i < ceil(N/16.0)*2; i++) {
		cmult256(x256[i], y256[i], &res_re, &res_im); 
		z256[i] = cpack256(res_re, res_im); 
	}
}

static inline void componentwise_multiply_real_scalar(int16_t *x,int16_t *y,int16_t *z, int N) {
	int i;
	for (i = 1; i <= N; i++) {
		z[i] = x[i] * y[i];
	}
}

static inline void componentwise_multiply_real_128(int16_t *x,int16_t *y,int16_t *z, int N) {

	__m128i *x128 = (__m128i *)x;
	__m128i *y128 = (__m128i *)y;
	__m128i *z128 = (__m128i *)z;

	int i;
	for (i = 0; i < ceil(N/8.0); i++) {
		z128[i] = _mm_mulhrs_epi16(x128[i], y128[i]); 
	}
}

static inline void componentwise_multiply_real_256(int16_t *x,int16_t *y,int16_t *z, uint16_t N) {

	__m256i *x256 = (__m256i *)x;
	__m256i *y256 = (__m256i *)y;
	__m256i *z256 = (__m256i *)z;

	int i;
	for (i = 0; i < ceil(N/16.0); i+=1) {
		z256[i] = _mm256_mulhrs_epi16(x256[i], y256[i]);
	}

}
