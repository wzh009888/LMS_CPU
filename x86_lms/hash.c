#include <string.h>
#include "hash.h"
#include "sha256.h"
#include "sha256_eight.h"
#include "sha256_sixteen.h"
#include "hss_zeroize.h"

#define ALLOW_VERBOSE 0         /* 1 -> we allow the dumping of intermediate */
                                /*      states.  Useful for debugging; horrid */
                                /*      for security */

/*
 * This is the file that implements the hashing APIs we use internally.
 * At the present, our parameter sets support only one hash function
 * (SHA-256, using full 256 bit output), however, that is likely to change
 * in the future
 */

#if ALLOW_VERBOSE
#include <stdio.h>
#include <stdbool.h>
/*
 * Debugging flag; if this is set, we chat about what we're hashing, and what
 * the result is it's useful when debugging; however we probably don't want to
 * do this if we're multithreaded...
 */
bool hss_verbose = false;
#endif

/*
 * This will hash the message, given the hash type. It assumes that the result
 * buffer is large enough for the hash
 */
void hss_hash_ctx(void *result, int hash_type, union hash_context *ctx,
		  const void *message, size_t message_len)
{
#if ALLOW_VERBOSE
	if (hss_verbose) {
		int i; for (i = 0; i < message_len; i++) printf(" %02x%s", ((unsigned char *)message)[i], (i % 16 == 15) ? "\n" : "");
	}
#endif

	switch (hash_type) {
	case HASH_SHA256: {
#if USE_OPENSSL == 1
		SHA256_Init(&ctx->sha256);
		SHA256_Update(&ctx->sha256, message, message_len);
		SHA256_Final(result, &ctx->sha256);
#else
		self_SHA256_Init(&ctx->sha256);
		self_SHA256_Update(&ctx->sha256, message, message_len);
		self_SHA256_Final(result, &ctx->sha256);
#endif
#if ALLOW_VERBOSE
		if (hss_verbose) {
			printf(" ->");
			int i; for (i = 0; i < 32; i++) printf(" %02x", ((unsigned char *)result)[i]); printf("\n");
		}
#endif
		break;
	}
	}
}

#if USE_OPENSSL != 1
//new version -- eight
void hss_hash_ctx_eight(void *result0, void *result1, void *result2,
			void *result3, void *result4, void *result5,
			void *result6, void *result7, int hash_type,
			const void *message0, const void *message1, const void *message2, const void *message3,
			const void *message4, const void *message5, const void *message6, const void *message7, size_t message_len)
{
	union hash_context ctx[8];

#if ALLOW_VERBOSE
	if (hss_verbose) {
		int i; for (i = 0; i < message_len; i++) printf(" %02x%s", ((unsigned char *)message)[i], (i % 16 == 15) ? "\n" : "");
	}
#endif

	switch (hash_type) {
	case HASH_SHA256: {
		self_SHA256_Init(&ctx[0].sha256);
		self_SHA256_Init(&ctx[1].sha256);
		self_SHA256_Init(&ctx[2].sha256);
		self_SHA256_Init(&ctx[3].sha256);
		self_SHA256_Init(&ctx[4].sha256);
		self_SHA256_Init(&ctx[5].sha256);
		self_SHA256_Init(&ctx[6].sha256);
		self_SHA256_Init(&ctx[7].sha256);

		eight_SHA256_Update(&ctx[0].sha256, &ctx[1].sha256, &ctx[2].sha256, &ctx[3].sha256, &ctx[4].sha256, &ctx[5].sha256, &ctx[6].sha256, &ctx[7].sha256, message0, message1, message2, message3, message4, message5, message6, message7, message_len);
		eight_SHA256_Final(result0, result1, result2, result3, result4, result5, result6, result7, &ctx[0].sha256, &ctx[1].sha256, &ctx[2].sha256, &ctx[3].sha256, &ctx[4].sha256, &ctx[5].sha256, &ctx[6].sha256, &ctx[7].sha256);
#if ALLOW_VERBOSE
		if (hss_verbose) {
			printf(" ->");
			int i; for (i = 0; i < 32; i++) printf(" %02x", ((unsigned char *)result)[i]); printf("\n");
		}
#endif
		break;
	}
	}
}
#endif
//end new version -- eight
#if USE_OPENSSL != 1
//new version -- sixteen
void hss_hash_ctx_sixteen(void *result0, void *result1, void *result2, void *result3,
			  void *result4, void *result5, void *result6, void *result7,
			  void *result8, void *result9, void *result10, void *result11,
			  void *result12, void *result13, void *result14, void *result15, int hash_type,
			  const void *message0, const void *message1, const void *message2, const void *message3,
			  const void *message4, const void *message5, const void *message6, const void *message7,
			  const void *message8, const void *message9, const void *message10, const void *message11,
			  const void *message12, const void *message13, const void *message14, const void *message15, size_t message_len)
{
	union hash_context ctx[16];

#if ALLOW_VERBOSE
	if (hss_verbose) {
		int i; for (i = 0; i < message_len; i++) printf(" %02x%s", ((unsigned char *)message)[i], (i % 16 == 15) ? "\n" : "");
	}
#endif

	switch (hash_type) {
	case HASH_SHA256: {
		self_SHA256_Init(&ctx[0].sha256);
		self_SHA256_Init(&ctx[1].sha256);
		self_SHA256_Init(&ctx[2].sha256);
		self_SHA256_Init(&ctx[3].sha256);
		self_SHA256_Init(&ctx[4].sha256);
		self_SHA256_Init(&ctx[5].sha256);
		self_SHA256_Init(&ctx[6].sha256);
		self_SHA256_Init(&ctx[7].sha256);
		self_SHA256_Init(&ctx[8].sha256);
		self_SHA256_Init(&ctx[9].sha256);
		self_SHA256_Init(&ctx[10].sha256);
		self_SHA256_Init(&ctx[11].sha256);
		self_SHA256_Init(&ctx[12].sha256);
		self_SHA256_Init(&ctx[13].sha256);
		self_SHA256_Init(&ctx[14].sha256);
		self_SHA256_Init(&ctx[15].sha256);
		sixteen_SHA256_Update(&ctx[0].sha256, &ctx[1].sha256, &ctx[2].sha256, &ctx[3].sha256, &ctx[4].sha256, &ctx[5].sha256, &ctx[6].sha256, &ctx[7].sha256, &ctx[8].sha256, &ctx[9].sha256, &ctx[10].sha256, &ctx[11].sha256, &ctx[12].sha256, &ctx[13].sha256, &ctx[14].sha256, &ctx[15].sha256, message0, message1, message2, message3, message4, message5, message6, message7, message8, message9, message10, message11, message12, message13, message14, message15, message_len);
		sixteen_SHA256_Final(result0, result1, result2, result3, result4, result5, result6, result7, result8, result9, result10, result11, result12, result13, result14, result15, &ctx[0].sha256, &ctx[1].sha256, &ctx[2].sha256, &ctx[3].sha256, &ctx[4].sha256, &ctx[5].sha256, &ctx[6].sha256, &ctx[7].sha256, &ctx[8].sha256, &ctx[9].sha256, &ctx[10].sha256, &ctx[11].sha256, &ctx[12].sha256, &ctx[13].sha256, &ctx[14].sha256, &ctx[15].sha256);
#if ALLOW_VERBOSE
		if (hss_verbose) {
			printf(" ->");
			int i; for (i = 0; i < 32; i++) printf(" %02x", ((unsigned char *)result)[i]); printf("\n");
		}
#endif
		break;
	}
	}
}
#endif
//end new version -- sixteen

void hss_hash(void *result, int hash_type,
	      const void *message, size_t message_len)
{
	union hash_context ctx;

	hss_hash_ctx(result, hash_type, &ctx, message, message_len);
	hss_zeroize(&ctx, sizeof ctx);
}


/*
 * This provides an API to do incremental hashing.  We use it when hashing the
 * message; since we don't know how long it could be, we don't want to
 * allocate a buffer that's long enough for that, plus the decoration we add
 */
void hss_init_hash_context(int h, union hash_context *ctx)
{
	switch (h) {
	case HASH_SHA256:
#if USE_OPENSSL == 1
		SHA256_Init(&ctx->sha256);
#else
		self_SHA256_Init(&ctx->sha256);
#endif
		break;
	}
}

void hss_update_hash_context(int h, union hash_context *ctx,
			     const void *msg, size_t len_msg)
{
#if ALLOW_VERBOSE
	if (hss_verbose) {
		int i; for (i = 0; i < len_msg; i++) printf(" %02x", ((unsigned char *)msg)[i]);
	}
#endif
	switch (h) {
	case HASH_SHA256:
#if USE_OPENSSL == 1
		SHA256_Update(&ctx->sha256, msg, len_msg);
#else
		self_SHA256_Update(&ctx->sha256, msg, len_msg);
#endif
		break;
	}
}

void hss_finalize_hash_context(int h, union hash_context *ctx, void *buffer)
{
	switch (h) {
	case HASH_SHA256:
#if USE_OPENSSL == 1
		SHA256_Final(buffer, &ctx->sha256);
#else
		self_SHA256_Final(buffer, &ctx->sha256);
#endif
#if ALLOW_VERBOSE
		if (hss_verbose) {
			printf(" -->");
			int i; for (i = 0; i < 32; i++) printf(" %02x", ((unsigned char *)buffer)[i]);
			printf("\n");
		}
#endif
		break;
	}
}


unsigned hss_hash_length(int hash_type)
{
	switch (hash_type) {
	case HASH_SHA256: return 32;
	}
	return 0;
}

unsigned hss_hash_blocksize(int hash_type)
{
	switch (hash_type) {
	case HASH_SHA256: return 64;
	}
	return 0;
}
