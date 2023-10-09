#include "sha256_sixteen.h"
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

# define Sigma0(x)  _mm512_xor_epi32(_mm512_xor_si512(_mm512_rol_epi32(x, 30), _mm512_rol_epi32(x, 19)), _mm512_rol_epi32(x, 10))

# define Sigma1(x)  _mm512_xor_si512(_mm512_xor_si512(_mm512_rol_epi32(x, 26), _mm512_rol_epi32(x, 21)), _mm512_rol_epi32(x, 7))

# define sigma0(x)  _mm512_xor_si512(_mm512_xor_si512(_mm512_rol_epi32(x, 25), _mm512_rol_epi32(x, 14)), _mm512_srli_epi32(x, 3))

# define sigma1(x)  _mm512_xor_si512(_mm512_xor_si512(_mm512_rol_epi32(x, 15), _mm512_rol_epi32(x, 13)), _mm512_srli_epi32(x, 10))

# define Ch(x, y, z)     (_mm512_xor_si512(z, _mm512_and_si512(x, _mm512_xor_si512(y, z))))
# define Maj(x, y, z)    _mm512_or_si512(_mm512_and_si512(y, _mm512_or_si512(x, z)), _mm512_and_si512(x, z))



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
__m512i kv16[64];

// __m256i byte_swap(__m256i tt_data)
// {
//     __m256i map, mdst;
//
//     map = _mm256_set_epi8(
// 	12, 13, 14, 15, 8, 9, 10, 11, 4, 5, 6, 7, 0, 1, 2, 3,   // first 128-bit lane
// 	12, 13, 14, 15, 8, 9, 10, 11, 4, 5, 6, 7, 0, 1, 2, 3 // second 128-bit lane
// 	);
//     mdst = _mm256_shuffle_epi8(tt_data, map);
// 	return mdst;
// }

#define FASTER
#ifdef FASTER

#define HOST_c21(i) do { \
		tt_data = _mm512_setr_epi32( \
			data[0][count * 16 + i], data[1][count * 16 + i], \
			data[2][count * 16 + i], data[3][count * 16 + i], \
			data[4][count * 16 + i], data[5][count * 16 + i], \
			data[6][count * 16 + i], data[7][count * 16 + i], \
			data[8][count * 16 + i], data[9][count * 16 + i], \
			data[10][count * 16 + i], data[11][count * 16 + i], \
			data[12][count * 16 + i], data[13][count * 16 + i], \
			data[14][count * 16 + i], data[15][count * 16 + i]); \
		T1v16 = lv16 = (_mm512_srli_epi32(tt_data, 24) \
				| _mm512_slli_epi32(tt_data, 24) \
				| _mm512_srli_epi32(_mm512_and_si512(tt_data, f1), 8) \
				| _mm512_slli_epi32(_mm512_and_si512(tt_data, f0), 8)); \
		Xv16[i] = lv16; } while(0)

#define ROUND_00_15(i, av16, bv16, cv16, dv16, ev16, fv16, gv16, hv16)          do {    \
		T1v16 = _mm512_add_epi32(T1v16, hv16); \
		T1v16 = _mm512_add_epi32(T1v16, Sigma1(ev16)); \
		T1v16 = _mm512_add_epi32(T1v16, Ch(ev16, fv16, gv16)); \
		T1v16 = _mm512_add_epi32(T1v16, kv16[i]); \
		hv16 = _mm512_add_epi32(Sigma0(av16), Maj(av16, bv16, cv16)); \
		dv16 = _mm512_add_epi32(dv16, T1v16); \
		hv16 = _mm512_add_epi32(hv16, T1v16); \
} while (0)

#define ROUND_16_63(i, a, b, c, d, e, f, g, h, X) \
	do {    \
		s0v16 = Xv16[(i + 1) & 0x0f];    s0v16 = sigma0(s0v16);        \
		s1v16 = Xv16[(i + 14) & 0x0f];   s1v16 = sigma1(s1v16);      \
		s0v16 = _mm512_add_epi32(s0v16, s1v16); \
		s0v16 = _mm512_add_epi32(Xv16[(i) & 0x0f], s0v16);      \
		Xv16[(i) & 0x0f] = _mm512_add_epi32(s0v16, Xv16[(i + 9) & 0x0f]);  \
		T1v16 = Xv16[(i) & 0x0f]; \
		ROUND_00_15(i, a, b, c, d, e, f, g, h); } while (0)

void sha256_block_data_order_sixteen(
	SHA256_CTX *ctx0, SHA256_CTX *ctx1, SHA256_CTX *ctx2,
	SHA256_CTX *ctx3, SHA256_CTX *ctx4, SHA256_CTX *ctx5,
	SHA256_CTX *ctx6, SHA256_CTX *ctx7, SHA256_CTX *ctx8,
	SHA256_CTX *ctx9, SHA256_CTX *ctx10, SHA256_CTX *ctx11,
	SHA256_CTX *ctx12, SHA256_CTX *ctx13, SHA256_CTX *ctx14, SHA256_CTX *ctx15,
	const void *in0, const void *in1, const void *in2, const void *in3,
	const void *in4, const void *in5, const void *in6, const void *in7,
	const void *in8, const void *in9, const void *in10, const void *in11,
	const void *in12, const void *in13, const void *in14, const void *in15,
	size_t num)
{
	SHA256_CTX *ctx[16];
	const void *in[16];
	unsigned int a16[16], b16[16], c16[16], d16[16], e16[16], f16[16], g16[16], h16[16],
		     s016[16], s116[16], T116[16], T216[16];
	__m512i av16, bv16, cv16, dv16, ev16, fv16, gv16, hv16, s0v16, s1v16, T1v16, T2v16;
	__m512i t_av16, t_bv16, t_cv16, t_dv16, t_ev16, t_fv16, t_gv16, t_hv16;
	unsigned int X16[16][16], l16[16];
	__m512i Xv16[16], lv16;
	int i, j;
	__m512i kv16[64];
	__m512i f0, f1;
	unsigned int fff0, fff1;
	const int *data[16];
	int *data_big[16];
	__m512i tt_data;

	for (i = 0; i < 64; i++) {
		kv16[i] = _mm512_set_epi32(
			K256[i], K256[i], K256[i], K256[i],
			K256[i], K256[i], K256[i], K256[i],
			K256[i], K256[i], K256[i], K256[i],
			K256[i], K256[i], K256[i], K256[i]);
	}
	ctx[0] = ctx0; ctx[1] = ctx1; ctx[2] = ctx2; ctx[3] = ctx3;
	ctx[4] = ctx4; ctx[5] = ctx5; ctx[6] = ctx6; ctx[7] = ctx7;
	ctx[8] = ctx8; ctx[9] = ctx9; ctx[10] = ctx10; ctx[11] = ctx11;
	ctx[12] = ctx12; ctx[13] = ctx13; ctx[14] = ctx14; ctx[15] = ctx15;
	in[0] = in0; in[1] = in1; in[2] = in2; in[3] = in3;
	in[4] = in4; in[5] = in5; in[6] = in6; in[7] = in7;
	in[8] = in8; in[9] = in9; in[10] = in10; in[11] = in11;
	in[12] = in12; in[13] = in13; in[14] = in14; in[15] = in15;


	fff0 = 0xff00;
	fff1 = 0xff0000;
	f0 = _mm512_set_epi32(
		fff0, fff0, fff0, fff0,
		fff0, fff0, fff0, fff0,
		fff0, fff0, fff0, fff0,
		fff0, fff0, fff0, fff0);
	f1 = _mm512_set_epi32(
		fff1, fff1, fff1, fff1,
		fff1, fff1, fff1, fff1,
		fff1, fff1, fff1, fff1,
		fff1, fff1, fff1, fff1);

	for (i = 0; i < 16; i++)
		data[i] = (const int *)in[i];
	int count = 0;

	av16 = _mm512_set_epi32(
		ctx[0]->h[0], ctx[1]->h[0], ctx[2]->h[0], ctx[3]->h[0],
		ctx[4]->h[0], ctx[5]->h[0], ctx[6]->h[0], ctx[7]->h[0],
		ctx[8]->h[0], ctx[9]->h[0], ctx[10]->h[0], ctx[11]->h[0],
		ctx[12]->h[0], ctx[13]->h[0], ctx[14]->h[0], ctx[15]->h[0]);
	bv16 = _mm512_set_epi32(
		ctx[0]->h[1], ctx[1]->h[1], ctx[2]->h[1], ctx[3]->h[1],
		ctx[4]->h[1], ctx[5]->h[1], ctx[6]->h[1], ctx[7]->h[1],
		ctx[8]->h[1], ctx[9]->h[1], ctx[10]->h[1], ctx[11]->h[1],
		ctx[12]->h[1], ctx[13]->h[1], ctx[14]->h[1], ctx[15]->h[1]);
	cv16 = _mm512_set_epi32(
		ctx[0]->h[2], ctx[1]->h[2], ctx[2]->h[2], ctx[3]->h[2],
		ctx[4]->h[2], ctx[5]->h[2], ctx[6]->h[2], ctx[7]->h[2],
		ctx[8]->h[2], ctx[9]->h[2], ctx[10]->h[2], ctx[11]->h[2],
		ctx[12]->h[2], ctx[13]->h[2], ctx[14]->h[2], ctx[15]->h[2]);
	dv16 = _mm512_set_epi32(
		ctx[0]->h[3], ctx[1]->h[3], ctx[2]->h[3], ctx[3]->h[3],
		ctx[4]->h[3], ctx[5]->h[3], ctx[6]->h[3], ctx[7]->h[3],
		ctx[8]->h[3], ctx[9]->h[3], ctx[10]->h[3], ctx[11]->h[3],
		ctx[12]->h[3], ctx[13]->h[3], ctx[14]->h[3], ctx[15]->h[3]);
	ev16 = _mm512_set_epi32(
		ctx[0]->h[4], ctx[1]->h[4], ctx[2]->h[4], ctx[3]->h[4],
		ctx[4]->h[4], ctx[5]->h[4], ctx[6]->h[4], ctx[7]->h[4],
		ctx[8]->h[4], ctx[9]->h[4], ctx[10]->h[4], ctx[11]->h[4],
		ctx[12]->h[4], ctx[13]->h[4], ctx[14]->h[4], ctx[15]->h[4]);
	fv16 = _mm512_set_epi32(
		ctx[0]->h[5], ctx[1]->h[5], ctx[2]->h[5], ctx[3]->h[5],
		ctx[4]->h[5], ctx[5]->h[5], ctx[6]->h[5], ctx[7]->h[5],
		ctx[8]->h[5], ctx[9]->h[5], ctx[10]->h[5], ctx[11]->h[5],
		ctx[12]->h[5], ctx[13]->h[5], ctx[14]->h[5], ctx[15]->h[5]);
	gv16 = _mm512_set_epi32(
		ctx[0]->h[6], ctx[1]->h[6], ctx[2]->h[6], ctx[3]->h[6],
		ctx[4]->h[6], ctx[5]->h[6], ctx[6]->h[6], ctx[7]->h[6],
		ctx[8]->h[6], ctx[9]->h[6], ctx[10]->h[6], ctx[11]->h[6],
		ctx[12]->h[6], ctx[13]->h[6], ctx[14]->h[6], ctx[15]->h[6]);
	hv16 = _mm512_set_epi32(
		ctx[0]->h[7], ctx[1]->h[7], ctx[2]->h[7], ctx[3]->h[7],
		ctx[4]->h[7], ctx[5]->h[7], ctx[6]->h[7], ctx[7]->h[7],
		ctx[8]->h[7], ctx[9]->h[7], ctx[10]->h[7], ctx[11]->h[7],
		ctx[12]->h[7], ctx[13]->h[7], ctx[14]->h[7], ctx[15]->h[7]);


	while (num--) {
		t_av16 = av16;
		t_bv16 = bv16;
		t_cv16 = cv16;
		t_dv16 = dv16;
		t_ev16 = ev16;
		t_fv16 = fv16;
		t_gv16 = gv16;
		t_hv16 = hv16;

		HOST_c21(0);
		ROUND_00_15(0, av16, bv16, cv16, dv16, ev16, fv16, gv16, hv16);
		HOST_c21(1);
		ROUND_00_15(1, hv16, av16, bv16, cv16, dv16, ev16, fv16, gv16);
		HOST_c21(2);
		ROUND_00_15(2, gv16, hv16, av16, bv16, cv16, dv16, ev16, fv16);
		HOST_c21(3);
		ROUND_00_15(3, fv16, gv16, hv16, av16, bv16, cv16, dv16, ev16);
		HOST_c21(4);
		ROUND_00_15(4, ev16, fv16, gv16, hv16, av16, bv16, cv16, dv16);
		HOST_c21(5);
		ROUND_00_15(5, dv16, ev16, fv16, gv16, hv16, av16, bv16, cv16);
		HOST_c21(6);
		ROUND_00_15(6, cv16, dv16, ev16, fv16, gv16, hv16, av16, bv16);
		HOST_c21(7);
		ROUND_00_15(7, bv16, cv16, dv16, ev16, fv16, gv16, hv16, av16);
		HOST_c21(8);
		ROUND_00_15(8, av16, bv16, cv16, dv16, ev16, fv16, gv16, hv16);
		HOST_c21(9);
		ROUND_00_15(9, hv16, av16, bv16, cv16, dv16, ev16, fv16, gv16);
		HOST_c21(10);
		ROUND_00_15(10, gv16, hv16, av16, bv16, cv16, dv16, ev16, fv16);
		HOST_c21(11);
		ROUND_00_15(11, fv16, gv16, hv16, av16, bv16, cv16, dv16, ev16);
		HOST_c21(12);
		ROUND_00_15(12, ev16, fv16, gv16, hv16, av16, bv16, cv16, dv16);
		HOST_c21(13);
		ROUND_00_15(13, dv16, ev16, fv16, gv16, hv16, av16, bv16, cv16);
		HOST_c21(14);
		ROUND_00_15(14, cv16, dv16, ev16, fv16, gv16, hv16, av16, bv16);
		HOST_c21(15);
		ROUND_00_15(15, bv16, cv16, dv16, ev16, fv16, gv16, hv16, av16);

		for (i = 16; i < 64; i += 8) {
			ROUND_16_63(i + 0, av16, bv16, cv16, dv16, ev16, fv16, gv16, hv16, Xv16);
			ROUND_16_63(i + 1, hv16, av16, bv16, cv16, dv16, ev16, fv16, gv16, Xv16);
			ROUND_16_63(i + 2, gv16, hv16, av16, bv16, cv16, dv16, ev16, fv16, Xv16);
			ROUND_16_63(i + 3, fv16, gv16, hv16, av16, bv16, cv16, dv16, ev16, Xv16);
			ROUND_16_63(i + 4, ev16, fv16, gv16, hv16, av16, bv16, cv16, dv16, Xv16);
			ROUND_16_63(i + 5, dv16, ev16, fv16, gv16, hv16, av16, bv16, cv16, Xv16);
			ROUND_16_63(i + 6, cv16, dv16, ev16, fv16, gv16, hv16, av16, bv16, Xv16);
			ROUND_16_63(i + 7, bv16, cv16, dv16, ev16, fv16, gv16, hv16, av16, Xv16);
		}

		av16 = _mm512_add_epi32(av16, t_av16);
		bv16 = _mm512_add_epi32(bv16, t_bv16);
		cv16 = _mm512_add_epi32(cv16, t_cv16);
		dv16 = _mm512_add_epi32(dv16, t_dv16);
		ev16 = _mm512_add_epi32(ev16, t_ev16);
		fv16 = _mm512_add_epi32(fv16, t_fv16);
		gv16 = _mm512_add_epi32(gv16, t_gv16);
		hv16 = _mm512_add_epi32(hv16, t_hv16);

		count++;
	}

	_mm512_store_si512((__m512i *)&a16[0], av16);
	_mm512_store_si512((__m512i *)&b16[0], bv16);
	_mm512_store_si512((__m512i *)&c16[0], cv16);
	_mm512_store_si512((__m512i *)&d16[0], dv16);
	_mm512_store_si512((__m512i *)&e16[0], ev16);
	_mm512_store_si512((__m512i *)&f16[0], fv16);
	_mm512_store_si512((__m512i *)&g16[0], gv16);
	_mm512_store_si512((__m512i *)&h16[0], hv16);

	for (i = 0; i < 16; i++) {
		ctx[i]->h[0] = a16[i];
		ctx[i]->h[1] = b16[i];
		ctx[i]->h[2] = c16[i];
		ctx[i]->h[3] = d16[i];
		ctx[i]->h[4] = e16[i];
		ctx[i]->h[5] = f16[i];
		ctx[i]->h[6] = g16[i];
		ctx[i]->h[7] = h16[i];
	}
}

#else
void sha256_block_data_order_sixteen(
	SHA256_CTX *ctx0, SHA256_CTX *ctx1, SHA256_CTX *ctx2,
	SHA256_CTX *ctx3, SHA256_CTX *ctx4, SHA256_CTX *ctx5,
	SHA256_CTX *ctx6, SHA256_CTX *ctx7, SHA256_CTX *ctx8,
	SHA256_CTX *ctx9, SHA256_CTX *ctx10, SHA256_CTX *ctx11,
	SHA256_CTX *ctx12, SHA256_CTX *ctx13, SHA256_CTX *ctx14, SHA256_CTX *ctx15,
	const void *in0, const void *in1, const void *in2, const void *in3,
	const void *in4, const void *in5, const void *in6, const void *in7,
	const void *in8, const void *in9, const void *in10, const void *in11,
	const void *in12, const void *in13, const void *in14, const void *in15,
	size_t num)
{
	SHA256_CTX *ctx[16];
	const void *in[16];
	unsigned int a16[16], b16[16], c16[16], d16[16], e16[16], f16[16], g16[16], h16[16],
		     s016[16], s116[16], T116[16], T216[16];
	__m512i av16, bv16, cv16, dv16, ev16, fv16, gv16, hv16, s0v16, s1v16, T1v16, T2v16;
	__m512i t_av16, t_bv16, t_cv16, t_dv16, t_ev16, t_fv16, t_gv16, t_hv16;
	unsigned int X16[16][16], l16[16];
	__m512i Xv16[16], lv16;
	int i, j;
	__m512i kv16[64];
	__m512i f0, f1;
	unsigned int fff0, fff1;
	const int *data[16];
	int *data_big[16];
	__m512i tt_data;

	for (i = 0; i < 64; i++) {
		kv16[i] = _mm512_set_epi32(
			K256[i], K256[i], K256[i], K256[i],
			K256[i], K256[i], K256[i], K256[i],
			K256[i], K256[i], K256[i], K256[i],
			K256[i], K256[i], K256[i], K256[i]);
	}
	ctx[0] = ctx0; ctx[1] = ctx1; ctx[2] = ctx2; ctx[3] = ctx3;
	ctx[4] = ctx4; ctx[5] = ctx5; ctx[6] = ctx6; ctx[7] = ctx7;
	ctx[8] = ctx8; ctx[9] = ctx9; ctx[10] = ctx10; ctx[11] = ctx11;
	ctx[12] = ctx12; ctx[13] = ctx13; ctx[14] = ctx14; ctx[15] = ctx15;
	in[0] = in0; in[1] = in1; in[2] = in2; in[3] = in3;
	in[4] = in4; in[5] = in5; in[6] = in6; in[7] = in7;
	in[8] = in8; in[9] = in9; in[10] = in10; in[11] = in11;
	in[12] = in12; in[13] = in13; in[14] = in14; in[15] = in15;

	fff0 = 0xff00;
	fff1 = 0xff0000;
	f0 = _mm512_set_epi32(
		fff0, fff0, fff0, fff0,
		fff0, fff0, fff0, fff0,
		fff0, fff0, fff0, fff0,
		fff0, fff0, fff0, fff0);
	f1 = _mm512_set_epi32(
		fff1, fff1, fff1, fff1,
		fff1, fff1, fff1, fff1,
		fff1, fff1, fff1, fff1,
		fff1, fff1, fff1, fff1);

	for (i = 0; i < 16; i++)
		data[i] = (const int *)in[i];
	int count = 0;

	av16 = _mm512_set_epi32(
		ctx[0]->h[0], ctx[1]->h[0], ctx[2]->h[0], ctx[3]->h[0],
		ctx[4]->h[0], ctx[5]->h[0], ctx[6]->h[0], ctx[7]->h[0],
		ctx[8]->h[0], ctx[9]->h[0], ctx[10]->h[0], ctx[11]->h[0],
		ctx[12]->h[0], ctx[13]->h[0], ctx[14]->h[0], ctx[15]->h[0]);
	bv16 = _mm512_set_epi32(
		ctx[0]->h[1], ctx[1]->h[1], ctx[2]->h[1], ctx[3]->h[1],
		ctx[4]->h[1], ctx[5]->h[1], ctx[6]->h[1], ctx[7]->h[1],
		ctx[8]->h[1], ctx[9]->h[1], ctx[10]->h[1], ctx[11]->h[1],
		ctx[12]->h[1], ctx[13]->h[1], ctx[14]->h[1], ctx[15]->h[1]);
	cv16 = _mm512_set_epi32(
		ctx[0]->h[2], ctx[1]->h[2], ctx[2]->h[2], ctx[3]->h[2],
		ctx[4]->h[2], ctx[5]->h[2], ctx[6]->h[2], ctx[7]->h[2],
		ctx[8]->h[2], ctx[9]->h[2], ctx[10]->h[2], ctx[11]->h[2],
		ctx[12]->h[2], ctx[13]->h[2], ctx[14]->h[2], ctx[15]->h[2]);
	dv16 = _mm512_set_epi32(
		ctx[0]->h[3], ctx[1]->h[3], ctx[2]->h[3], ctx[3]->h[3],
		ctx[4]->h[3], ctx[5]->h[3], ctx[6]->h[3], ctx[7]->h[3],
		ctx[8]->h[3], ctx[9]->h[3], ctx[10]->h[3], ctx[11]->h[3],
		ctx[12]->h[3], ctx[13]->h[3], ctx[14]->h[3], ctx[15]->h[3]);
	ev16 = _mm512_set_epi32(
		ctx[0]->h[4], ctx[1]->h[4], ctx[2]->h[4], ctx[3]->h[4],
		ctx[4]->h[4], ctx[5]->h[4], ctx[6]->h[4], ctx[7]->h[4],
		ctx[8]->h[4], ctx[9]->h[4], ctx[10]->h[4], ctx[11]->h[4],
		ctx[12]->h[4], ctx[13]->h[4], ctx[14]->h[4], ctx[15]->h[4]);
	fv16 = _mm512_set_epi32(
		ctx[0]->h[5], ctx[1]->h[5], ctx[2]->h[5], ctx[3]->h[5],
		ctx[4]->h[5], ctx[5]->h[5], ctx[6]->h[5], ctx[7]->h[5],
		ctx[8]->h[5], ctx[9]->h[5], ctx[10]->h[5], ctx[11]->h[5],
		ctx[12]->h[5], ctx[13]->h[5], ctx[14]->h[5], ctx[15]->h[5]);
	gv16 = _mm512_set_epi32(
		ctx[0]->h[6], ctx[1]->h[6], ctx[2]->h[6], ctx[3]->h[6],
		ctx[4]->h[6], ctx[5]->h[6], ctx[6]->h[6], ctx[7]->h[6],
		ctx[8]->h[6], ctx[9]->h[6], ctx[10]->h[6], ctx[11]->h[6],
		ctx[12]->h[6], ctx[13]->h[6], ctx[14]->h[6], ctx[15]->h[6]);
	hv16 = _mm512_set_epi32(
		ctx[0]->h[7], ctx[1]->h[7], ctx[2]->h[7], ctx[3]->h[7],
		ctx[4]->h[7], ctx[5]->h[7], ctx[6]->h[7], ctx[7]->h[7],
		ctx[8]->h[7], ctx[9]->h[7], ctx[10]->h[7], ctx[11]->h[7],
		ctx[12]->h[7], ctx[13]->h[7], ctx[14]->h[7], ctx[15]->h[7]);

	while (num--) {
		t_av16 = av16;
		t_bv16 = bv16;
		t_cv16 = cv16;
		t_dv16 = dv16;
		t_ev16 = ev16;
		t_fv16 = fv16;
		t_gv16 = gv16;
		t_hv16 = hv16;

		for (i = 0; i < 16; i++) {
			tt_data = _mm512_setr_epi32(
				data[0][count * 16 + i], data[1][count * 16 + i],
				data[2][count * 16 + i], data[3][count * 16 + i],
				data[4][count * 16 + i], data[5][count * 16 + i],
				data[6][count * 16 + i], data[7][count * 16 + i],
				data[8][count * 16 + i], data[9][count * 16 + i],
				data[10][count * 16 + i], data[11][count * 16 + i],
				data[12][count * 16 + i], data[13][count * 16 + i],
				data[14][count * 16 + i], data[15][count * 16 + i]);

			lv16 = (_mm512_srli_epi32(tt_data, 24)
				| _mm512_slli_epi32(tt_data, 24)
				| _mm512_srli_epi32(_mm512_and_si512(tt_data, f1), 8)
				| _mm512_slli_epi32(_mm512_and_si512(tt_data, f0), 8));
			//lv8 = byte_swap(tt_data);

			Xv16[i] = lv16;
			T1v16 = _mm512_add_epi32(lv16, hv16);
			T1v16 = _mm512_add_epi32(T1v16, Sigma1(ev16));
			T1v16 = _mm512_add_epi32(T1v16, Ch(ev16, fv16, gv16));
			T1v16 = _mm512_add_epi32(T1v16, kv16[i]);

			T2v16 = _mm512_add_epi32(Sigma0(av16), Maj(av16, bv16, cv16));
			hv16 = gv16;
			gv16 = fv16;
			fv16 = ev16;
			ev16 = _mm512_add_epi32(dv16, T1v16);
			dv16 = cv16;
			cv16 = bv16;
			bv16 = av16;
			av16 = _mm512_add_epi32(T1v16, T2v16);
		}

		for (i = 16; i < 64; i++) {
			s0v16 = Xv16[(i + 1) & 0x0f];
			s0v16 = sigma0(s0v16);
			s1v16 = Xv16[(i + 14) & 0x0f];
			s1v16 = sigma1(s1v16);
			Xv16[i & 0xf] = _mm512_add_epi32(Xv16[i & 0xf], s0v16);
			Xv16[i & 0xf] = _mm512_add_epi32(Xv16[i & 0xf], s1v16);
			Xv16[i & 0xf] = _mm512_add_epi32(Xv16[i & 0xf], Xv16[(i + 9) & 0xf]);
			T1v16 = _mm512_add_epi32(Xv16[i & 0xf], hv16);
			T1v16 = _mm512_add_epi32(T1v16, Sigma1(ev16));
			T1v16 = _mm512_add_epi32(T1v16, Ch(ev16, fv16, gv16));
			T1v16 = _mm512_add_epi32(T1v16, kv16[i]);
			T2v16 = _mm512_add_epi32(Sigma0(av16), Maj(av16, bv16, cv16));
			hv16 = gv16;
			gv16 = fv16;
			fv16 = ev16;
			ev16 = _mm512_add_epi32(dv16, T1v16);
			dv16 = cv16;
			cv16 = bv16;
			bv16 = av16;
			av16 = _mm512_add_epi32(T1v16, T2v16);
		}

		av16 = _mm512_add_epi32(av16, t_av16);
		bv16 = _mm512_add_epi32(bv16, t_bv16);
		cv16 = _mm512_add_epi32(cv16, t_cv16);
		dv16 = _mm512_add_epi32(dv16, t_dv16);
		ev16 = _mm512_add_epi32(ev16, t_ev16);
		fv16 = _mm512_add_epi32(fv16, t_fv16);
		gv16 = _mm512_add_epi32(gv16, t_gv16);
		hv16 = _mm512_add_epi32(hv16, t_hv16);

		count++;
	}

	_mm512_store_si512((__m512i *)&a16[0], av16);
	_mm512_store_si512((__m512i *)&b16[0], bv16);
	_mm512_store_si512((__m512i *)&c16[0], cv16);
	_mm512_store_si512((__m512i *)&d16[0], dv16);
	_mm512_store_si512((__m512i *)&e16[0], ev16);
	_mm512_store_si512((__m512i *)&f16[0], fv16);
	_mm512_store_si512((__m512i *)&g16[0], gv16);
	_mm512_store_si512((__m512i *)&h16[0], hv16);

	for (i = 0; i < 16; i++) {
		ctx[i]->h[0] = a16[i];
		ctx[i]->h[1] = b16[i];
		ctx[i]->h[2] = c16[i];
		ctx[i]->h[3] = d16[i];
		ctx[i]->h[4] = e16[i];
		ctx[i]->h[5] = f16[i];
		ctx[i]->h[6] = g16[i];
		ctx[i]->h[7] = h16[i];
	}
}
#endif

int sixteen_SHA256_Update(
	SHA256_CTX *c0, SHA256_CTX *c1,
	SHA256_CTX *c2, SHA256_CTX *c3,
	SHA256_CTX *c4, SHA256_CTX *c5,
	SHA256_CTX *c6, SHA256_CTX *c7,
	SHA256_CTX *c8, SHA256_CTX *c9,
	SHA256_CTX *c10, SHA256_CTX *c11,
	SHA256_CTX *c12, SHA256_CTX *c13,
	SHA256_CTX *c14, SHA256_CTX *c15,
	const unsigned char *data0, const unsigned char *data1,
	const unsigned char *data2, const unsigned char *data3,
	const unsigned char *data4, const unsigned char *data5,
	const unsigned char *data6, const unsigned char *data7,
	const unsigned char *data8, const unsigned char *data9,
	const unsigned char *data10, const unsigned char *data11,
	const unsigned char *data12, const unsigned char *data13,
	const unsigned char *data14, const unsigned char *data15,
	size_t len)
{
	unsigned char *p;
	unsigned int l;
	size_t n;

	if (len == 0)
		return 1;

	l = (((unsigned int)len) << 3) & 0xffffffffUL;
	if (l < 0)                              /* overflow */
		c0->Nh = c1->Nh = c2->Nh = c3->Nh = c4->Nh =
			c5->Nh = c6->Nh = c7->Nh = c8->Nh = c9->Nh =
				c10->Nh = c11->Nh = c12->Nh = c13->Nh = c14->Nh = c15->Nh = 1;
	c0->Nh += (unsigned int)(len >> 29);
	c1->Nh += (unsigned int)(len >> 29);
	c2->Nh += (unsigned int)(len >> 29);
	c3->Nh += (unsigned int)(len >> 29);
	c4->Nh += (unsigned int)(len >> 29);
	c5->Nh += (unsigned int)(len >> 29);
	c6->Nh += (unsigned int)(len >> 29);
	c7->Nh += (unsigned int)(len >> 29);
	c8->Nh += (unsigned int)(len >> 29);
	c9->Nh += (unsigned int)(len >> 29);
	c10->Nh += (unsigned int)(len >> 29);
	c11->Nh += (unsigned int)(len >> 29);
	c12->Nh += (unsigned int)(len >> 29);
	c13->Nh += (unsigned int)(len >> 29);
	c14->Nh += (unsigned int)(len >> 29);
	c15->Nh += (unsigned int)(len >> 29);

	c0->Nl = c1->Nl = c2->Nl = c3->Nl =
		c4->Nl = c5->Nl = c6->Nl = c7->Nl =
			c8->Nl = c9->Nl = c10->Nl = c11->Nl =
				c12->Nl = c13->Nl = c14->Nl = c15->Nl = l;

	n = len / HASH_CBLOCK;
	if (n > 0) {
		sha256_block_data_order_sixteen(
			c0, c1, c2, c3, c4, c5, c6, c7,
			c8, c9, c10, c11, c12, c13, c14, c15,
			data0, data1, data2, data3, data4, data5, data6, data7,
			data8, data9, data10, data11, data12, data13, data14, data15, n);
		// printf("data0 = %02x\n", data0);
		n *= HASH_CBLOCK;
		data0 += n;
		data1 += n;
		data2 += n;
		data3 += n;
		data4 += n;
		data5 += n;
		data6 += n;
		data7 += n;
		data8 += n;
		data9 += n;
		data10 += n;
		data11 += n;
		data12 += n;
		data13 += n;
		data14 += n;
		data15 += n;
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
		p = (unsigned char *)c8->data;
		c8->num = (unsigned int)len;
		memcpy(p, data8, len);
		p = (unsigned char *)c9->data;
		c9->num = (unsigned int)len;
		memcpy(p, data9, len);
		p = (unsigned char *)c10->data;
		c10->num = (unsigned int)len;
		memcpy(p, data10, len);
		p = (unsigned char *)c11->data;
		c11->num = (unsigned int)len;
		memcpy(p, data11, len);
		p = (unsigned char *)c12->data;
		c12->num = (unsigned int)len;
		memcpy(p, data12, len);
		p = (unsigned char *)c13->data;
		c13->num = (unsigned int)len;
		memcpy(p, data13, len);
		p = (unsigned char *)c14->data;
		c14->num = (unsigned int)len;
		memcpy(p, data14, len);
		p = (unsigned char *)c15->data;
		c15->num = (unsigned int)len;
		memcpy(p, data15, len);
	}
	return 1;
}

#define HASH_MAKE_STRING(c, s)   do {    \
		unsigned long ll;               \
		unsigned int nn;               \
		for (nn = 0; nn < SHA256_DIGEST_LENGTH / 4; nn++)       \
		{ ll = (c)->h[nn]; (void)HOST_l2c(ll, (s)); }  \
} while (0)

int sixteen_SHA256_Final(
	unsigned char *md0, unsigned char *md1,
	unsigned char *md2, unsigned char *md3,
	unsigned char *md4, unsigned char *md5,
	unsigned char *md6, unsigned char *md7,
	unsigned char *md8, unsigned char *md9,
	unsigned char *md10, unsigned char *md11,
	unsigned char *md12, unsigned char *md13,
	unsigned char *md14, unsigned char *md15,
	SHA256_CTX *c0, SHA256_CTX *c1,
	SHA256_CTX *c2, SHA256_CTX *c3,
	SHA256_CTX *c4, SHA256_CTX *c5,
	SHA256_CTX *c6, SHA256_CTX *c7,
	SHA256_CTX *c8, SHA256_CTX *c9,
	SHA256_CTX *c10, SHA256_CTX *c11,
	SHA256_CTX *c12, SHA256_CTX *c13,
	SHA256_CTX *c14, SHA256_CTX *c15)
{
	unsigned char *p0 = (unsigned char *)c0->data;
	unsigned char *p1 = (unsigned char *)c1->data;
	unsigned char *p2 = (unsigned char *)c2->data;
	unsigned char *p3 = (unsigned char *)c3->data;
	unsigned char *p4 = (unsigned char *)c4->data;
	unsigned char *p5 = (unsigned char *)c5->data;
	unsigned char *p6 = (unsigned char *)c6->data;
	unsigned char *p7 = (unsigned char *)c7->data;
	unsigned char *p8 = (unsigned char *)c8->data;
	unsigned char *p9 = (unsigned char *)c9->data;
	unsigned char *p10 = (unsigned char *)c10->data;
	unsigned char *p11 = (unsigned char *)c11->data;
	unsigned char *p12 = (unsigned char *)c12->data;
	unsigned char *p13 = (unsigned char *)c13->data;
	unsigned char *p14 = (unsigned char *)c14->data;
	unsigned char *p15 = (unsigned char *)c15->data;
	size_t n = c1->num;

	p0[n] = p1[n] = p2[n] = p3[n] = p4[n] = p5[n] = p6[n] = p7[n] = p8[n] = p9[n] = p10[n] = p11[n] = p12[n] = p13[n] = p14[n] = p15[n] = 0x80;
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
		memset(p8 + n, 0, HASH_CBLOCK - n);
		memset(p9 + n, 0, HASH_CBLOCK - n);
		memset(p10 + n, 0, HASH_CBLOCK - n);
		memset(p11 + n, 0, HASH_CBLOCK - n);
		memset(p12 + n, 0, HASH_CBLOCK - n);
		memset(p13 + n, 0, HASH_CBLOCK - n);
		memset(p14 + n, 0, HASH_CBLOCK - n);
		memset(p15 + n, 0, HASH_CBLOCK - n);
		n = 0;
		sha256_block_data_order_sixteen(
			c0, c1, c2, c3, c4, c5, c6, c7,
			c8, c9, c10, c11, c12, c13, c14, c15,
			p0, p1, p2, p3, p4, p5, p6, p7,
			p8, p9, p10, p11, p12, p13, p14, p15, 1);
	}
	memset(p0 + n, 0, HASH_CBLOCK - 8 - n);
	memset(p1 + n, 0, HASH_CBLOCK - 8 - n);
	memset(p2 + n, 0, HASH_CBLOCK - 8 - n);
	memset(p3 + n, 0, HASH_CBLOCK - 8 - n);
	memset(p4 + n, 0, HASH_CBLOCK - 8 - n);
	memset(p5 + n, 0, HASH_CBLOCK - 8 - n);
	memset(p6 + n, 0, HASH_CBLOCK - 8 - n);
	memset(p7 + n, 0, HASH_CBLOCK - 8 - n);
	memset(p8 + n, 0, HASH_CBLOCK - 8 - n);
	memset(p9 + n, 0, HASH_CBLOCK - 8 - n);
	memset(p10 + n, 0, HASH_CBLOCK - 8 - n);
	memset(p11 + n, 0, HASH_CBLOCK - 8 - n);
	memset(p12 + n, 0, HASH_CBLOCK - 8 - n);
	memset(p13 + n, 0, HASH_CBLOCK - 8 - n);
	memset(p14 + n, 0, HASH_CBLOCK - 8 - n);
	memset(p15 + n, 0, HASH_CBLOCK - 8 - n);

	p0 += HASH_CBLOCK - 8;
	p1 += HASH_CBLOCK - 8;
	p2 += HASH_CBLOCK - 8;
	p3 += HASH_CBLOCK - 8;
	p4 += HASH_CBLOCK - 8;
	p5 += HASH_CBLOCK - 8;
	p6 += HASH_CBLOCK - 8;
	p7 += HASH_CBLOCK - 8;
	p8 += HASH_CBLOCK - 8;
	p9 += HASH_CBLOCK - 8;
	p10 += HASH_CBLOCK - 8;
	p11 += HASH_CBLOCK - 8;
	p12 += HASH_CBLOCK - 8;
	p13 += HASH_CBLOCK - 8;
	p14 += HASH_CBLOCK - 8;
	p15 += HASH_CBLOCK - 8;

	(void)HOST_l2c(c0->Nh, p0);
	(void)HOST_l2c(c1->Nh, p1);
	(void)HOST_l2c(c2->Nh, p2);
	(void)HOST_l2c(c3->Nh, p3);
	(void)HOST_l2c(c4->Nh, p4);
	(void)HOST_l2c(c5->Nh, p5);
	(void)HOST_l2c(c6->Nh, p6);
	(void)HOST_l2c(c7->Nh, p7);
	(void)HOST_l2c(c8->Nh, p8);
	(void)HOST_l2c(c9->Nh, p9);
	(void)HOST_l2c(c10->Nh, p10);
	(void)HOST_l2c(c11->Nh, p11);
	(void)HOST_l2c(c12->Nh, p12);
	(void)HOST_l2c(c13->Nh, p13);
	(void)HOST_l2c(c14->Nh, p14);
	(void)HOST_l2c(c15->Nh, p15);

	(void)HOST_l2c(c0->Nl, p0);
	(void)HOST_l2c(c1->Nl, p1);
	(void)HOST_l2c(c2->Nl, p2);
	(void)HOST_l2c(c3->Nl, p3);
	(void)HOST_l2c(c4->Nl, p4);
	(void)HOST_l2c(c5->Nl, p5);
	(void)HOST_l2c(c6->Nl, p6);
	(void)HOST_l2c(c7->Nl, p7);
	(void)HOST_l2c(c8->Nl, p8);
	(void)HOST_l2c(c9->Nl, p9);
	(void)HOST_l2c(c10->Nl, p10);
	(void)HOST_l2c(c11->Nl, p11);
	(void)HOST_l2c(c12->Nl, p12);
	(void)HOST_l2c(c13->Nl, p13);
	(void)HOST_l2c(c14->Nl, p14);
	(void)HOST_l2c(c15->Nl, p15);

	p0 -= HASH_CBLOCK;
	p1 -= HASH_CBLOCK;
	p2 -= HASH_CBLOCK;
	p3 -= HASH_CBLOCK;
	p4 -= HASH_CBLOCK;
	p5 -= HASH_CBLOCK;
	p6 -= HASH_CBLOCK;
	p7 -= HASH_CBLOCK;
	p8 -= HASH_CBLOCK;
	p9 -= HASH_CBLOCK;
	p10 -= HASH_CBLOCK;
	p11 -= HASH_CBLOCK;
	p12 -= HASH_CBLOCK;
	p13 -= HASH_CBLOCK;
	p14 -= HASH_CBLOCK;
	p15 -= HASH_CBLOCK;

	sha256_block_data_order_sixteen(
		c0, c1, c2, c3, c4, c5, c6, c7,
		c8, c9, c10, c11, c12, c13, c14, c15,
		p0, p1, p2, p3, p4, p5, p6, p7,
		p8, p9, p10, p11, p12, p13, p14, p15, 1);

	HASH_MAKE_STRING(c0, md0);
	HASH_MAKE_STRING(c1, md1);
	HASH_MAKE_STRING(c2, md2);
	HASH_MAKE_STRING(c3, md3);
	HASH_MAKE_STRING(c4, md4);
	HASH_MAKE_STRING(c5, md5);
	HASH_MAKE_STRING(c6, md6);
	HASH_MAKE_STRING(c7, md7);
	HASH_MAKE_STRING(c8, md8);
	HASH_MAKE_STRING(c9, md9);
	HASH_MAKE_STRING(c10, md10);
	HASH_MAKE_STRING(c11, md11);
	HASH_MAKE_STRING(c12, md12);
	HASH_MAKE_STRING(c13, md13);
	HASH_MAKE_STRING(c14, md14);
	HASH_MAKE_STRING(c15, md15);

	return 1;
}
