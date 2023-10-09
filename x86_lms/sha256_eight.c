#include "sha256_eight.h"
#include <immintrin.h>
#include <arpa/inet.h>

#define HASH_LONG               SHA_LONG
#define HASH_CBLOCK             SHA_CBLOCK
# define SHA_LONG               unsigned int
# define SHA_CBLOCK             (SHA_LBLOCK * 4)

#define self_HASH_UPDATE             self_SHA256_Update
#define self_HASH_FINAL              self_SHA256_Final
#define self_HASH_CTX                SHA256_CTX

# define HOST_l2c(l, c)  (*((c)++) = (unsigned char)(((l) >> 24) & 0xff),      \
			  *((c)++) = (unsigned char)(((l) >> 16) & 0xff),      \
			  *((c)++) = (unsigned char)(((l) >> 8) & 0xff),      \
			  *((c)++) = (unsigned char)(((l)) & 0xff),      \
			  l)

# define Sigma0(x)  _mm256_xor_si256(_mm256_xor_si256(_mm256_rol_epi32(x, 30), _mm256_rol_epi32(x, 19)), _mm256_rol_epi32(x, 10))

# define Sigma1(x)  _mm256_xor_si256(_mm256_xor_si256(_mm256_rol_epi32(x, 26), _mm256_rol_epi32(x, 21)), _mm256_rol_epi32(x, 7))

# define sigma0(x)  _mm256_xor_si256(_mm256_xor_si256(_mm256_rol_epi32(x, 25), _mm256_rol_epi32(x, 14)), _mm256_srli_epi32(x, 3))

# define sigma1(x)  _mm256_xor_si256(_mm256_xor_si256(_mm256_rol_epi32(x, 15), _mm256_rol_epi32(x, 13)), _mm256_srli_epi32(x, 10))

# define Ch(x, y, z)     (_mm256_xor_si256(z,_mm256_and_si256(x,_mm256_xor_si256(y,z))))
# define Maj(x, y, z)    _mm256_or_si256(_mm256_and_si256(y, _mm256_or_si256(x,z)), _mm256_and_si256(x,z))


// # define ROTATE(v, n) (((v) << (n)) | ((v) >> (32 - (n))))
//
// # define Sigma0(x)       (ROTATE((x), 30) ^ ROTATE((x), 19) ^ ROTATE((x), 10))
// # define Sigma1(x)       (ROTATE((x), 26) ^ ROTATE((x), 21) ^ ROTATE((x), 7))
// # define sigma0(x)       (ROTATE((x), 25) ^ ROTATE((x), 14) ^ ((x) >> 3))
// # define sigma1(x)       (ROTATE((x), 15) ^ ROTATE((x), 13) ^ ((x) >> 10))
//
// # define Ch(x, y, z)     ((z) ^ ((x) & ((y) ^ (z))))
// # define Maj(x, y, z)    (((y) & ((x) | (z))) | (x) & (z))


static SHA_LONG K256[64] = {
	0x428a2f98UL, 0x71374491UL, 0xb5c0fbcfUL, 0xe9b5dba5UL,
	0x3956c25bUL, 0x59f111f1UL, 0x923f82a4UL, 0xab1c5ed5UL,
	0xd807aa98UL, 0x12835b01UL, 0x243185beUL, 0x550c7dc3UL,
	0x72be5d74UL, 0x80deb1feUL, 0x9bdc06a7UL, 0xc19bf174UL,
	0xe49b69c1UL, 0xefbe4786UL, 0x0fc19dc6UL, 0x240ca1ccUL,
	0x2de92c6fUL, 0x4a7484aaUL, 0x5cb0a9dcUL, 0x76f988daUL,
	0x983e5152UL, 0xa831c66dUL, 0xb00327c8UL, 0xbf597fc7UL,
	0xc6e00bf3UL, 0xd5a79147UL, 0x06ca6351UL, 0x14292967UL,
	0x27b70a85UL, 0x2e1b2138UL, 0x4d2c6dfcUL, 0x53380d13UL,
	0x650a7354UL, 0x766a0abbUL, 0x81c2c92eUL, 0x92722c85UL,
	0xa2bfe8a1UL, 0xa81a664bUL, 0xc24b8b70UL, 0xc76c51a3UL,
	0xd192e819UL, 0xd6990624UL, 0xf40e3585UL, 0x106aa070UL,
	0x19a4c116UL, 0x1e376c08UL, 0x2748774cUL, 0x34b0bcb5UL,
	0x391c0cb3UL, 0x4ed8aa4aUL, 0x5b9cca4fUL, 0x682e6ff3UL,
	0x748f82eeUL, 0x78a5636fUL, 0x84c87814UL, 0x8cc70208UL,
	0x90befffaUL, 0xa4506cebUL, 0xbef9a3f7UL, 0xc67178f2UL
};
__m256i kv8[64];

__m256i byte_swap(__m256i tt_data)
{
    __m256i map, mdst;

    map = _mm256_set_epi8(
	12, 13, 14, 15, 8, 9, 10, 11, 4, 5, 6, 7, 0, 1, 2, 3,   // first 128-bit lane
	12, 13, 14, 15, 8, 9, 10, 11, 4, 5, 6, 7, 0, 1, 2, 3 // second 128-bit lane
	);
    mdst = _mm256_shuffle_epi8(tt_data, map);
	return mdst;
}

#define FASTER
#ifdef FASTER

#define HOST_c21(i) do { \
	tt_data = _mm256_setr_epi32( \
			data[0][count * 16 + i], data[1][count * 16 + i], \
			data[2][count * 16 + i], data[3][count * 16 + i], \
			data[4][count * 16 + i], data[5][count * 16 + i], \
			data[6][count * 16 + i], data[7][count * 16 + i]); \
	T1v8 = lv8 = byte_swap(tt_data); \
	Xv8[i] = lv8; } while(0)

#define ROUND_00_15(i, av8, bv8, cv8, dv8, ev8, fv8, gv8, hv8)          do {    \
		T1v8 = _mm256_add_epi32(T1v8, hv8); \
		T1v8 = _mm256_add_epi32(T1v8, Sigma1(ev8)); \
		T1v8 = _mm256_add_epi32(T1v8, Ch(ev8, fv8, gv8)); \
		T1v8 = _mm256_add_epi32(T1v8, kv8[i]); \
		hv8 = _mm256_add_epi32(Sigma0(av8), Maj(av8, bv8, cv8)); \
		dv8 = _mm256_add_epi32(dv8, T1v8); \
		hv8 = _mm256_add_epi32(hv8, T1v8); } while (0)

#define ROUND_16_63(i, a, b, c, d, e, f, g, h, X) \
		do {    \
		s0v8 = Xv8[(i + 1) & 0x0f];    s0v8 = sigma0(s0v8);        \
		s1v8 = Xv8[(i + 14) & 0x0f];   s1v8 = sigma1(s1v8);      \
		s0v8 = _mm256_add_epi32(s0v8, s1v8);	\
		s0v8 = _mm256_add_epi32(Xv8[(i) & 0x0f], s0v8);	\
		Xv8[(i) & 0x0f] = _mm256_add_epi32(s0v8, Xv8[(i + 9) & 0x0f]);  \
		T1v8 = Xv8[(i) & 0x0f]; \
		ROUND_00_15(i, a, b, c, d, e, f, g, h); } while (0)

void sha256_block_data_order_eight(
	SHA256_CTX *ctx0, SHA256_CTX *ctx1, SHA256_CTX *ctx2,
	SHA256_CTX *ctx3, SHA256_CTX *ctx4, SHA256_CTX *ctx5,
	SHA256_CTX *ctx6, SHA256_CTX *ctx7,
	const void *in0, const void *in1, const void *in2, const void *in3,
	const void *in4, const void *in5, const void *in6, const void *in7,
	size_t num)
{
	SHA256_CTX *ctx[8];
	const void *in[8];
	unsigned int a8[8], b8[8], c8[8], d8[8], e8[8], f8[8], g8[8], h8[8],
		     s08[8], s18[8], T18[8], T28[8];
	__m256i av8, bv8, cv8, dv8, ev8, fv8, gv8, hv8, s0v8, s1v8, T1v8, T2v8;
	__m256i t_av8, t_bv8, t_cv8, t_dv8, t_ev8, t_fv8, t_gv8, t_hv8;
	unsigned int X8[16][8], l8[8];
	__m256i Xv8[16], lv8;
	int i, j;

	__m256i f0, f1;
	unsigned int fff0, fff1;
	const int *data[8];
	int *data_big[8];
	__m256i tt_data;

	for (i = 0; i < 64; i++) {
		kv8[i] = _mm256_set_epi32(
			K256[i], K256[i], K256[i], K256[i],
			K256[i], K256[i], K256[i], K256[i]);
	}
	ctx[0] = ctx0; ctx[1] = ctx1; ctx[2] = ctx2; ctx[3] = ctx3;
	ctx[4] = ctx4; ctx[5] = ctx5; ctx[6] = ctx6; ctx[7] = ctx7;
	in[0] = in0; in[1] = in1; in[2] = in2; in[3] = in3;
	in[4] = in4; in[5] = in5; in[6] = in6; in[7] = in7;


	fff0 = 0xff00;
	fff1 = 0xff0000;
	f0 = _mm256_set_epi32(
		fff0, fff0, fff0, fff0,
		fff0, fff0, fff0, fff0);
	f1 = _mm256_set_epi32(
		fff1, fff1, fff1, fff1,
		fff1, fff1, fff1, fff1);

	for (i = 0; i < 8; i++){
		data[i] = (const int *)in[i];
	}
	int count = 0;

	av8 = _mm256_set_epi32(
		ctx[0]->h[0], ctx[1]->h[0], ctx[2]->h[0], ctx[3]->h[0],
		ctx[4]->h[0], ctx[5]->h[0], ctx[6]->h[0], ctx[7]->h[0]);
	bv8 = _mm256_set_epi32(
		ctx[0]->h[1], ctx[1]->h[1], ctx[2]->h[1], ctx[3]->h[1],
		ctx[4]->h[1], ctx[5]->h[1], ctx[6]->h[1], ctx[7]->h[1]);
	cv8 = _mm256_set_epi32(
		ctx[0]->h[2], ctx[1]->h[2], ctx[2]->h[2], ctx[3]->h[2],
		ctx[4]->h[2], ctx[5]->h[2], ctx[6]->h[2], ctx[7]->h[2]);
	dv8 = _mm256_set_epi32(
		ctx[0]->h[3], ctx[1]->h[3], ctx[2]->h[3], ctx[3]->h[3],
		ctx[4]->h[3], ctx[5]->h[3], ctx[6]->h[3], ctx[7]->h[3]);
	ev8 = _mm256_set_epi32(
		ctx[0]->h[4], ctx[1]->h[4], ctx[2]->h[4], ctx[3]->h[4],
		ctx[4]->h[4], ctx[5]->h[4], ctx[6]->h[4], ctx[7]->h[4]);
	fv8 = _mm256_set_epi32(
		ctx[0]->h[5], ctx[1]->h[5], ctx[2]->h[5], ctx[3]->h[5],
		ctx[4]->h[5], ctx[5]->h[5], ctx[6]->h[5], ctx[7]->h[5]);
	gv8 = _mm256_set_epi32(
		ctx[0]->h[6], ctx[1]->h[6], ctx[2]->h[6], ctx[3]->h[6],
		ctx[4]->h[6], ctx[5]->h[6], ctx[6]->h[6], ctx[7]->h[6]);
	hv8 = _mm256_set_epi32(
		ctx[0]->h[7], ctx[1]->h[7], ctx[2]->h[7], ctx[3]->h[7],
		ctx[4]->h[7], ctx[5]->h[7], ctx[6]->h[7], ctx[7]->h[7]);

	while (num--) {
		t_av8 = av8;
		t_bv8 = bv8;
		t_cv8 = cv8;
		t_dv8 = dv8;
		t_ev8 = ev8;
		t_fv8 = fv8;
		t_gv8 = gv8;
		t_hv8 = hv8;

		HOST_c21(0);
		ROUND_00_15(0, av8, bv8, cv8, dv8, ev8, fv8, gv8, hv8);
		HOST_c21(1);
		ROUND_00_15(1, hv8, av8, bv8, cv8, dv8, ev8, fv8, gv8);
		HOST_c21(2);
		ROUND_00_15(2, gv8, hv8, av8, bv8, cv8, dv8, ev8, fv8);
		HOST_c21(3);
		ROUND_00_15(3, fv8, gv8, hv8, av8, bv8, cv8, dv8, ev8);
		HOST_c21(4);
		ROUND_00_15(4, ev8, fv8, gv8, hv8, av8, bv8, cv8, dv8);
		HOST_c21(5);
		ROUND_00_15(5, dv8, ev8, fv8, gv8, hv8, av8, bv8, cv8);
		HOST_c21(6);
		ROUND_00_15(6, cv8, dv8, ev8, fv8, gv8, hv8, av8, bv8);
		HOST_c21(7);
		ROUND_00_15(7, bv8, cv8, dv8, ev8, fv8, gv8, hv8, av8);
		HOST_c21(8);
		ROUND_00_15(8, av8, bv8, cv8, dv8, ev8, fv8, gv8, hv8);
		HOST_c21(9);
		ROUND_00_15(9, hv8, av8, bv8, cv8, dv8, ev8, fv8, gv8);
		HOST_c21(10);
		ROUND_00_15(10, gv8, hv8, av8, bv8, cv8, dv8, ev8, fv8);
		HOST_c21(11);
		ROUND_00_15(11, fv8, gv8, hv8, av8, bv8, cv8, dv8, ev8);
		HOST_c21(12);
		ROUND_00_15(12, ev8, fv8, gv8, hv8, av8, bv8, cv8, dv8);
		HOST_c21(13);
		ROUND_00_15(13, dv8, ev8, fv8, gv8, hv8, av8, bv8, cv8);
		HOST_c21(14);
		ROUND_00_15(14, cv8, dv8, ev8, fv8, gv8, hv8, av8, bv8);
		HOST_c21(15);
		ROUND_00_15(15, bv8, cv8, dv8, ev8, fv8, gv8, hv8, av8);

		for (i = 16; i < 64; i += 8) {
			ROUND_16_63(i + 0, av8, bv8, cv8, dv8, ev8, fv8, gv8, hv8, Xv8);
			ROUND_16_63(i + 1, hv8, av8, bv8, cv8, dv8, ev8, fv8, gv8, Xv8);
			ROUND_16_63(i + 2, gv8, hv8, av8, bv8, cv8, dv8, ev8, fv8, Xv8);
			ROUND_16_63(i + 3, fv8, gv8, hv8, av8, bv8, cv8, dv8, ev8, Xv8);
			ROUND_16_63(i + 4, ev8, fv8, gv8, hv8, av8, bv8, cv8, dv8, Xv8);
			ROUND_16_63(i + 5, dv8, ev8, fv8, gv8, hv8, av8, bv8, cv8, Xv8);
			ROUND_16_63(i + 6, cv8, dv8, ev8, fv8, gv8, hv8, av8, bv8, Xv8);
			ROUND_16_63(i + 7, bv8, cv8, dv8, ev8, fv8, gv8, hv8, av8, Xv8);
		}

		av8 = _mm256_add_epi32(av8, t_av8);
		bv8 = _mm256_add_epi32(bv8, t_bv8);
		cv8 = _mm256_add_epi32(cv8, t_cv8);
		dv8 = _mm256_add_epi32(dv8, t_dv8);
		ev8 = _mm256_add_epi32(ev8, t_ev8);
		fv8 = _mm256_add_epi32(fv8, t_fv8);
		gv8 = _mm256_add_epi32(gv8, t_gv8);
		hv8 = _mm256_add_epi32(hv8, t_hv8);

		count++;
	}

	_mm256_store_si256((__m256i *)&a8[0], av8);
	_mm256_store_si256((__m256i *)&b8[0], bv8);
	_mm256_store_si256((__m256i *)&c8[0], cv8);
	_mm256_store_si256((__m256i *)&d8[0], dv8);
	_mm256_store_si256((__m256i *)&e8[0], ev8);
	_mm256_store_si256((__m256i *)&f8[0], fv8);
	_mm256_store_si256((__m256i *)&g8[0], gv8);
	_mm256_store_si256((__m256i *)&h8[0], hv8);

	for (i = 0; i < 8; i++) {
		ctx[i]->h[0] = a8[i];
		ctx[i]->h[1] = b8[i];
		ctx[i]->h[2] = c8[i];
		ctx[i]->h[3] = d8[i];
		ctx[i]->h[4] = e8[i];
		ctx[i]->h[5] = f8[i];
		ctx[i]->h[6] = g8[i];
		ctx[i]->h[7] = h8[i];
	}
	// _mm256_storeu_si256((__m256i *)&s08[0], av8);
	// printf("%02x %02x %02x %02x %02x %02x %02x %02x\n", s08[0], s08[1], s08[2], s08[3], s08[4], s08[5], s08[6], s08[7]);
	//
	// printf("%02x %02x\n", a8[0], a8[1]);
	// //
	// int aa = 0;
	// while(1){
	// 	aa++;
	// }
	// printf("aa = %d\n", aa);

}

#else
void sha256_block_data_order_eight(
	SHA256_CTX *ctx0, SHA256_CTX *ctx1, SHA256_CTX *ctx2,
	SHA256_CTX *ctx3, SHA256_CTX *ctx4, SHA256_CTX *ctx5,
	SHA256_CTX *ctx6, SHA256_CTX *ctx7,
	const void *in0, const void *in1, const void *in2, const void *in3,
	const void *in4, const void *in5, const void *in6, const void *in7,
	size_t num)
{
	SHA256_CTX *ctx[8];
	const void *in[8];
	unsigned int a8[8], b8[8], c8[8], d8[8], e8[8], f8[8], g8[8], h8[8],
		     s08[8], s18[8], T18[8], T28[8];
	__m256i av8, bv8, cv8, dv8, ev8, fv8, gv8, hv8, s0v8, s1v8, T1v8, T2v8;
	__m256i t_av8, t_bv8, t_cv8, t_dv8, t_ev8, t_fv8, t_gv8, t_hv8;
	unsigned int X8[16][8], l8[8];
	__m256i Xv8[16], lv8;
	int i, j;
	__m256i kv8[64];
	__m256i f0, f1;
	unsigned int fff0, fff1;
	const int *data[8];
	int *data_big[8];
	__m256i tt_data;

	for (i = 0; i < 64; i++) {
		kv8[i] = _mm256_set_epi32(
			K256[i], K256[i], K256[i], K256[i],
			K256[i], K256[i], K256[i], K256[i]);
	}
	ctx[0] = ctx0; ctx[1] = ctx1; ctx[2] = ctx2; ctx[3] = ctx3;
	ctx[4] = ctx4; ctx[5] = ctx5; ctx[6] = ctx6; ctx[7] = ctx7;
	in[0] = in0; in[1] = in1; in[2] = in2; in[3] = in3;
	in[4] = in4; in[5] = in5; in[6] = in6; in[7] = in7;


	fff0 = 0xff00;
	fff1 = 0xff0000;
	f0 = _mm256_set_epi32(
		fff0, fff0, fff0, fff0,
		fff0, fff0, fff0, fff0);
	f1 = _mm256_set_epi32(
		fff1, fff1, fff1, fff1,
		fff1, fff1, fff1, fff1);

	for (i = 0; i < 8; i++){
		data[i] = (const int *)in[i];
	}
	int count = 0;

	av8 = _mm256_set_epi32(
		ctx[0]->h[0], ctx[1]->h[0], ctx[2]->h[0], ctx[3]->h[0],
		ctx[4]->h[0], ctx[5]->h[0], ctx[6]->h[0], ctx[7]->h[0]);
	bv8 = _mm256_set_epi32(
		ctx[0]->h[1], ctx[1]->h[1], ctx[2]->h[1], ctx[3]->h[1],
		ctx[4]->h[1], ctx[5]->h[1], ctx[6]->h[1], ctx[7]->h[1]);
	cv8 = _mm256_set_epi32(
		ctx[0]->h[2], ctx[1]->h[2], ctx[2]->h[2], ctx[3]->h[2],
		ctx[4]->h[2], ctx[5]->h[2], ctx[6]->h[2], ctx[7]->h[2]);
	dv8 = _mm256_set_epi32(
		ctx[0]->h[3], ctx[1]->h[3], ctx[2]->h[3], ctx[3]->h[3],
		ctx[4]->h[3], ctx[5]->h[3], ctx[6]->h[3], ctx[7]->h[3]);
	ev8 = _mm256_set_epi32(
		ctx[0]->h[4], ctx[1]->h[4], ctx[2]->h[4], ctx[3]->h[4],
		ctx[4]->h[4], ctx[5]->h[4], ctx[6]->h[4], ctx[7]->h[4]);
	fv8 = _mm256_set_epi32(
		ctx[0]->h[5], ctx[1]->h[5], ctx[2]->h[5], ctx[3]->h[5],
		ctx[4]->h[5], ctx[5]->h[5], ctx[6]->h[5], ctx[7]->h[5]);
	gv8 = _mm256_set_epi32(
		ctx[0]->h[6], ctx[1]->h[6], ctx[2]->h[6], ctx[3]->h[6],
		ctx[4]->h[6], ctx[5]->h[6], ctx[6]->h[6], ctx[7]->h[6]);
	hv8 = _mm256_set_epi32(
		ctx[0]->h[7], ctx[1]->h[7], ctx[2]->h[7], ctx[3]->h[7],
		ctx[4]->h[7], ctx[5]->h[7], ctx[6]->h[7], ctx[7]->h[7]);

	while (num--) {
		t_av8 = av8;
		t_bv8 = bv8;
		t_cv8 = cv8;
		t_dv8 = dv8;
		t_ev8 = ev8;
		t_fv8 = fv8;
		t_gv8 = gv8;
		t_hv8 = hv8;

		for (i = 0; i < 16; i++) {
			tt_data = _mm256_setr_epi32(
				data[0][count * 16 + i], data[1][count * 16 + i],
				data[2][count * 16 + i], data[3][count * 16 + i],
				data[4][count * 16 + i], data[5][count * 16 + i],
				data[6][count * 16 + i], data[7][count * 16 + i]);

			// lv8 = ( _mm256_srli_epi32(tt_data, 24)
			// 		| _mm256_slli_epi32(tt_data, 24)
			//       	| _mm256_srli_epi32(_mm256_and_si256(tt_data, f1), 8)
			// 		| _mm256_slli_epi32(_mm256_and_si256(tt_data, f0), 8));
			lv8 = byte_swap(tt_data);

			Xv8[i] = lv8;
			T1v8 = _mm256_add_epi32(lv8, hv8);
			T1v8 = _mm256_add_epi32(T1v8, Sigma1(ev8));
			T1v8 = _mm256_add_epi32(T1v8, Ch(ev8, fv8, gv8));
			T1v8 = _mm256_add_epi32(T1v8, kv8[i]);
			T2v8 = _mm256_add_epi32(Sigma0(av8), Maj(av8, bv8, cv8));
			hv8 = gv8;
			gv8 = fv8;
			fv8 = ev8;
			ev8 = _mm256_add_epi32(dv8, T1v8);
			dv8 = cv8;
			cv8 = bv8;
			bv8 = av8;
			av8 = _mm256_add_epi32(T1v8, T2v8);
			//
			// _mm256_storeu_si256((__m256i *)&s08[0], tt_data);
			// printf("%d tt_data %02x %02x %02x %02x %02x %02x %02x %02x\n", i, s08[0], s08[1], s08[2], s08[3], s08[4], s08[5], s08[6], s08[7]);

		}


		for (i = 16; i < 64; i++) {
			s0v8 = Xv8[(i + 1) & 0x0f];
			s0v8 = sigma0(s0v8);
			s1v8 = Xv8[(i + 14) & 0x0f];
			s1v8 = sigma1(s1v8);
			Xv8[i & 0xf] = _mm256_add_epi32(Xv8[i & 0xf], s0v8);
			Xv8[i & 0xf] = _mm256_add_epi32(Xv8[i & 0xf], s1v8);
			Xv8[i & 0xf] = _mm256_add_epi32(Xv8[i & 0xf], Xv8[(i + 9) & 0xf]);
			T1v8 = _mm256_add_epi32(Xv8[i & 0xf], hv8);
			T1v8 = _mm256_add_epi32(T1v8, Sigma1(ev8));
			T1v8 = _mm256_add_epi32(T1v8, Ch(ev8, fv8, gv8));
			T1v8 = _mm256_add_epi32(T1v8, kv8[i]);
			T2v8 = _mm256_add_epi32(Sigma0(av8), Maj(av8, bv8, cv8));
			hv8 = gv8;
			gv8 = fv8;
			fv8 = ev8;
			ev8 = _mm256_add_epi32(dv8, T1v8);
			dv8 = cv8;
			cv8 = bv8;
			bv8 = av8;
			av8 = _mm256_add_epi32(T1v8, T2v8);

			// _mm256_storeu_si256((__m256i *)&s08[0], av8);
			// printf("%d av8 %02x %02x %02x %02x %02x %02x %02x %02x\n", i, s08[0], s08[1], s08[2], s08[3], s08[4], s08[5], s08[6], s08[7]);
		}

		av8 = _mm256_add_epi32(av8, t_av8);
		bv8 = _mm256_add_epi32(bv8, t_bv8);
		cv8 = _mm256_add_epi32(cv8, t_cv8);
		dv8 = _mm256_add_epi32(dv8, t_dv8);
		ev8 = _mm256_add_epi32(ev8, t_ev8);
		fv8 = _mm256_add_epi32(fv8, t_fv8);
		gv8 = _mm256_add_epi32(gv8, t_gv8);
		hv8 = _mm256_add_epi32(hv8, t_hv8);

		count++;
	}

	_mm256_store_si256((__m256i *)&a8[0], av8);
	_mm256_store_si256((__m256i *)&b8[0], bv8);
	_mm256_store_si256((__m256i *)&c8[0], cv8);
	_mm256_store_si256((__m256i *)&d8[0], dv8);
	_mm256_store_si256((__m256i *)&e8[0], ev8);
	_mm256_store_si256((__m256i *)&f8[0], fv8);
	_mm256_store_si256((__m256i *)&g8[0], gv8);
	_mm256_store_si256((__m256i *)&h8[0], hv8);

	for (i = 0; i < 8; i++) {
		ctx[i]->h[0] = a8[i];
		ctx[i]->h[1] = b8[i];
		ctx[i]->h[2] = c8[i];
		ctx[i]->h[3] = d8[i];
		ctx[i]->h[4] = e8[i];
		ctx[i]->h[5] = f8[i];
		ctx[i]->h[6] = g8[i];
		ctx[i]->h[7] = h8[i];
	}


}
#endif

int eight_SHA256_Update(
	SHA256_CTX *c0, SHA256_CTX *c1,
	SHA256_CTX *c2, SHA256_CTX *c3,
	SHA256_CTX *c4, SHA256_CTX *c5,
	SHA256_CTX *c6, SHA256_CTX *c7,
	const unsigned char *data0, const unsigned char *data1,
	const unsigned char *data2, const unsigned char *data3,
	const unsigned char *data4, const unsigned char *data5,
	const unsigned char *data6, const unsigned char *data7,
	size_t len)
{
	unsigned char *p;
	unsigned int l;
	size_t n;

	if (len == 0)
		return 1;

	l = (((unsigned int)len) << 3) & 0xffffffffUL;
	if (l < 0)                              /* overflow */
		c0->Nh = c1->Nh = c2->Nh = c3->Nh =
			c4->Nh = c5->Nh = c6->Nh = c7->Nh = 1;
	c0->Nh += (unsigned int)(len >> 29);
	c1->Nh += (unsigned int)(len >> 29);
	c2->Nh += (unsigned int)(len >> 29);
	c3->Nh += (unsigned int)(len >> 29);
	c4->Nh += (unsigned int)(len >> 29);
	c5->Nh += (unsigned int)(len >> 29);
	c6->Nh += (unsigned int)(len >> 29);
	c7->Nh += (unsigned int)(len >> 29);
	c0->Nl = c1->Nl = c2->Nl = c3->Nl =
		c4->Nl = c5->Nl = c6->Nl = c7->Nl = l;

	n = len / HASH_CBLOCK;
	if (n > 0) {
		sha256_block_data_order_eight(
			c0, c1, c2, c3, c4, c5, c6, c7,
			data0, data1, data2, data3, data4, data5, data6, data7, n);

		n *= HASH_CBLOCK;
		data0 += n;
		data1 += n;
		data2 += n;
		data3 += n;
		data4 += n;
		data5 += n;
		data6 += n;
		data7 += n;
		len -= n;
	}

	if (len != 0) {
		p = (unsigned char *)c0->data;
		c0->num = (unsigned int)len;
		memcpy(p, data0, len);
		p = (unsigned char *)c1->data;
		c1->num = (unsigned int)len;
		memcpy(p, data1, len);
		p = (unsigned char *)c2->data;
		c2->num = (unsigned int)len;
		memcpy(p, data2, len);
		p = (unsigned char *)c3->data;
		c3->num = (unsigned int)len;
		memcpy(p, data3, len);
		p = (unsigned char *)c4->data;
		c4->num = (unsigned int)len;
		memcpy(p, data4, len);
		p = (unsigned char *)c5->data;
		c5->num = (unsigned int)len;
		memcpy(p, data5, len);
		p = (unsigned char *)c6->data;
		c6->num = (unsigned int)len;
		memcpy(p, data6, len);
		p = (unsigned char *)c7->data;
		c7->num = (unsigned int)len;
		memcpy(p, data7, len);
	}
	return 1;
}

#define HASH_MAKE_STRING(c, s)   do {    \
		unsigned long ll;               \
		unsigned int nn;               \
		for (nn = 0; nn < SHA256_DIGEST_LENGTH / 4; nn++)       \
		{ ll = (c)->h[nn]; (void)HOST_l2c(ll, (s)); }  \
} while (0)

int eight_SHA256_Final(
	unsigned char *md0, unsigned char *md1,
	unsigned char *md2, unsigned char *md3,
	unsigned char *md4, unsigned char *md5,
	unsigned char *md6, unsigned char *md7,
	SHA256_CTX *c0, SHA256_CTX *c1,
	SHA256_CTX *c2, SHA256_CTX *c3,
	SHA256_CTX *c4, SHA256_CTX *c5,
	SHA256_CTX *c6, SHA256_CTX *c7)
{
	unsigned char *p0 = (unsigned char *)c0->data;
	unsigned char *p1 = (unsigned char *)c1->data;
	unsigned char *p2 = (unsigned char *)c2->data;
	unsigned char *p3 = (unsigned char *)c3->data;
	unsigned char *p4 = (unsigned char *)c4->data;
	unsigned char *p5 = (unsigned char *)c5->data;
	unsigned char *p6 = (unsigned char *)c6->data;
	unsigned char *p7 = (unsigned char *)c7->data;
	size_t n = c1->num;

	p0[n] = p1[n] = p2[n] = p3[n] = p4[n] = p5[n] = p6[n] = p7[n] = 0x80;
	n++;

	if (n > (HASH_CBLOCK - 8)) {
		memset(p0 + n, 0, HASH_CBLOCK - n);
		memset(p1 + n, 0, HASH_CBLOCK - n);
		memset(p2 + n, 0, HASH_CBLOCK - n);
		memset(p3 + n, 0, HASH_CBLOCK - n);
		memset(p4 + n, 0, HASH_CBLOCK - n);
		memset(p5 + n, 0, HASH_CBLOCK - n);
		memset(p6 + n, 0, HASH_CBLOCK - n);
		memset(p7 + n, 0, HASH_CBLOCK - n);
		n = 0;
		sha256_block_data_order_eight(
			c0, c1, c2, c3, c4, c5, c6, c7,
			p0, p1, p2, p3, p4, p5, p6, p7, 1);
	}
	memset(p0 + n, 0, HASH_CBLOCK - 8 - n);
	memset(p1 + n, 0, HASH_CBLOCK - 8 - n);
	memset(p2 + n, 0, HASH_CBLOCK - 8 - n);
	memset(p3 + n, 0, HASH_CBLOCK - 8 - n);
	memset(p4 + n, 0, HASH_CBLOCK - 8 - n);
	memset(p5 + n, 0, HASH_CBLOCK - 8 - n);
	memset(p6 + n, 0, HASH_CBLOCK - 8 - n);
	memset(p7 + n, 0, HASH_CBLOCK - 8 - n);

	p0 += HASH_CBLOCK - 8;
	p1 += HASH_CBLOCK - 8;
	p2 += HASH_CBLOCK - 8;
	p3 += HASH_CBLOCK - 8;
	p4 += HASH_CBLOCK - 8;
	p5 += HASH_CBLOCK - 8;
	p6 += HASH_CBLOCK - 8;
	p7 += HASH_CBLOCK - 8;

	(void)HOST_l2c(c0->Nh, p0);
	(void)HOST_l2c(c1->Nh, p1);
	(void)HOST_l2c(c2->Nh, p2);
	(void)HOST_l2c(c3->Nh, p3);
	(void)HOST_l2c(c4->Nh, p4);
	(void)HOST_l2c(c5->Nh, p5);
	(void)HOST_l2c(c6->Nh, p6);
	(void)HOST_l2c(c7->Nh, p7);

	(void)HOST_l2c(c0->Nl, p0);
	(void)HOST_l2c(c1->Nl, p1);
	(void)HOST_l2c(c2->Nl, p2);
	(void)HOST_l2c(c3->Nl, p3);
	(void)HOST_l2c(c4->Nl, p4);
	(void)HOST_l2c(c5->Nl, p5);
	(void)HOST_l2c(c6->Nl, p6);
	(void)HOST_l2c(c7->Nl, p7);

	p0 -= HASH_CBLOCK;
	p1 -= HASH_CBLOCK;
	p2 -= HASH_CBLOCK;
	p3 -= HASH_CBLOCK;
	p4 -= HASH_CBLOCK;
	p5 -= HASH_CBLOCK;
	p6 -= HASH_CBLOCK;
	p7 -= HASH_CBLOCK;

	sha256_block_data_order_eight(
		c0, c1, c2, c3, c4, c5, c6, c7,
		p0, p1, p2, p3, p4, p5, p6, p7, 1);

	HASH_MAKE_STRING(c0, md0);
	HASH_MAKE_STRING(c1, md1);
	HASH_MAKE_STRING(c2, md2);
	HASH_MAKE_STRING(c3, md3);
	HASH_MAKE_STRING(c4, md4);
	HASH_MAKE_STRING(c5, md5);
	HASH_MAKE_STRING(c6, md6);
	HASH_MAKE_STRING(c7, md7);

	return 1;
}
