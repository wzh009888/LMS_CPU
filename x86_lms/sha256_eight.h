#include <string.h>
#include <stdio.h>
#include "sha256.h"


#define USING_OPT

# define SHA256_DIGEST_LENGTH    32
# define SHA_LBLOCK      16
# define SHA_LONG unsigned int
# define MD32_REG_T int

// typedef struct self_SHA256state_st {
// 	unsigned int	h[8];
// 	unsigned int	Nl, Nh;
// 	unsigned int	data[SHA_LBLOCK];
// 	unsigned int	num, md_len;
// } SHA256_CTX;

// int self_SHA256_Init(SHA256_CTX *c);

int eight_SHA256_Update(
	SHA256_CTX *c0, SHA256_CTX *c1,
	SHA256_CTX *c2, SHA256_CTX *c3,
	SHA256_CTX *c4, SHA256_CTX *c5,
	SHA256_CTX *c6, SHA256_CTX *c7,
	const unsigned char *data0, const unsigned char *data1,
	const unsigned char *data2, const unsigned char *data3,
	const unsigned char *data4, const unsigned char *data5,
	const unsigned char *data6, const unsigned char *data7,
	size_t len);

int eight_SHA256_Final(
	unsigned char *md0, unsigned char *md1,
	unsigned char *md2, unsigned char *md3,
	unsigned char *md4, unsigned char *md5,
	unsigned char *md6, unsigned char *md7,
	SHA256_CTX *c0, SHA256_CTX *c1,
	SHA256_CTX *c2, SHA256_CTX *c3,
	SHA256_CTX *c4, SHA256_CTX *c5,
	SHA256_CTX *c6, SHA256_CTX *c7);
