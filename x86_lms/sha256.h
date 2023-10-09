#if !defined(SHA256_H_)
#define SHA256_H_

#if defined(EXT_SHA256_H)
#include EXT_SHA256_H
#else

#include "para.h"

#if USE_OPENSSL

#include <openssl/sha.h>

#else

/* SHA256 context. */
// typedef struct {
//   unsigned long int h[8];            /* state; this is in the CPU native format */
//   unsigned long Nl, Nh;              /* number of bits processed so far */
//   unsigned num;                      /* number of bytes within the below */
//                                      /* buffer */
//   unsigned char data[64];            /* input buffer.  This is in byte vector format */
// } SHA256_CTX;
//
// void SHA256_Init(SHA256_CTX *);  /* context */
//
// void SHA256_Update(SHA256_CTX *, /* context */
//                   const void *, /* input block */
//                   unsigned int);/* length of input block */
//
// void SHA256_Final(unsigned char *,
//                  SHA256_CTX *);

//new version -- one
# define SHA256_DIGEST_LENGTH    32
# define SHA_LBLOCK      16
# define SHA_LONG unsigned int
# define MD32_REG_T int

typedef struct self_SHA256state_st {
	unsigned int	h[8];
	unsigned int	Nl, Nh;
	unsigned int	data[SHA_LBLOCK];
	unsigned int	num, md_len;
} SHA256_CTX;

int self_SHA256_Init(SHA256_CTX *c);

int self_SHA256_Update(SHA256_CTX *c, const void *data_, size_t len);

int self_SHA256_Final(unsigned char *md, SHA256_CTX *c);


#endif

#endif /* EXT_SHA256_H */

#if !defined(SHA256_LEN)
#define SHA256_LEN 32    /* The length of a SHA256 hash output */
#endif

#endif /* ifdef(SHA256_H_) */
