/*
 * This is the code that implements the one-time-signature part of the LMS hash
 * based signatures
 */
#include <string.h>
#include <stdio.h>
#include <time.h>
#include "common_defs.h"
#include "lm_ots.h"
#include "lm_ots_common.h"
#include "hash.h"
#include "sha256.h"
#include "endian.h"
#include "hss_zeroize.h"
#include "hss_derive.h"
#include "hss_internal.h"
#include "mpi.h"

double sum_time2 = 0;
double sum_time4 = 0;
extern int global_level_nodes;
extern MPI_Comm local_comm02;
extern MPI_Comm local_comm;
extern int num01;
#if (defined(vectorization))

#if (defined(vector_512))
/*vectorization -- sixteen */
bool mpi_lm_ots_generate_public_key(
	param_set_t lm_ots_type,
	const unsigned char *I, /* Public key identifier */
	merkle_index_t q,       /* Diversification string, 4 bytes value */
	struct seed_derive *seed,
	unsigned char *public_key, size_t public_key_len)
{
	int me, proc_num;

	MPI_Comm_rank(MPI_COMM_WORLD, &me);             //进程号
	MPI_Comm_size(MPI_COMM_WORLD, &proc_num);       //进程数
	// double start, end;

	// start = MPI_Wtime();
	unsigned int new_size;

	new_size = proc_num / global_level_nodes;

	// new_size = 2;
	MPI_Barrier(MPI_COMM_WORLD);

	// end = MPI_Wtime();
	// if(me == 0) printf("0 %lf\n", 1000 * (end - start));
	// start = MPI_Wtime();

	/* Look up the parameter set */
	unsigned h, n, w, p, ls;

	if (!lm_ots_look_up_parameter_set(lm_ots_type, &h, &n, &w, &p, &ls))
		return false;
	MPI_Barrier(MPI_COMM_WORLD);
	/* Start the hash that computes the final value */
	union hash_context public_ctx;

	hss_init_hash_context(h, &public_ctx);
	{
		unsigned char prehash_prefix[PBLC_PREFIX_LEN];
		memcpy(prehash_prefix + PBLC_I, I, I_LEN);
		put_bigendian(prehash_prefix + PBLC_Q, q, 4);
		SET_D(prehash_prefix + PBLC_D, D_PBLC);
		hss_update_hash_context(h, &public_ctx, prehash_prefix,
					PBLC_PREFIX_LEN);
	}
	int id;                                 //子通信域中新的进程号

#ifdef Gather_communicator
	MPI_Comm_rank(local_comm, &id);         //进程号
#else
	//假如没有划分通信子域
	id = me / global_level_nodes;
#endif
	// new_size是子通信域的总进程数
	int e = p / new_size;
	int f = p % new_size;
	unsigned int local = e + ((id < f) ? 1 : 0);
	unsigned int offset = id * e + ((id < f) ? id : f);

	hss_seed_derive_set_j(seed, 0);
	int num = local / 16;
	int rem = local % 16;
	int m_num = 16;
	unsigned char buffer[((p / 16) + 1) * 16][ITER_MAX_LEN];
	unsigned char *bufs[m_num];
	union hash_context ctx[16];
	unsigned char empty[ITER_MAX_LEN];

	for (size_t i = 0; i < m_num; i++)
		bufs[i] = &empty[0];

	for (size_t i = 0; i < p; i++) {
		memcpy(buffer[i] + ITER_I, I, I_LEN);
		put_bigendian(buffer[i] + ITER_Q, q, 4);
		hss_seed_derive(buffer[i] + ITER_PREV, seed, i < p - 1);
		put_bigendian(buffer[i] + ITER_K, i, 2);
	}
	for (size_t i = 0; i < num; i++) {
		int suboffset = offset + i * 16;
		for (size_t r = 0; r < m_num; r++)
			bufs[r] = buffer[suboffset + r];
		for (size_t j = 0; j < (1 << w) - 1; j++) {
			bufs[0][ITER_J] = j;
			bufs[1][ITER_J] = j;
			bufs[2][ITER_J] = j;
			bufs[3][ITER_J] = j;
			bufs[4][ITER_J] = j;
			bufs[5][ITER_J] = j;
			bufs[6][ITER_J] = j;
			bufs[7][ITER_J] = j;
			bufs[8][ITER_J] = j;
			bufs[9][ITER_J] = j;
			bufs[10][ITER_J] = j;
			bufs[11][ITER_J] = j;
			bufs[12][ITER_J] = j;
			bufs[13][ITER_J] = j;
			bufs[14][ITER_J] = j;
			bufs[15][ITER_J] = j;
			hss_hash_ctx_sixteen(
				bufs[0] + ITER_PREV, bufs[1] + ITER_PREV, bufs[2] + ITER_PREV, bufs[3] + ITER_PREV,
				bufs[4] + ITER_PREV, bufs[5] + ITER_PREV, bufs[6] + ITER_PREV, bufs[7] + ITER_PREV,
				bufs[8] + ITER_PREV, bufs[9] + ITER_PREV, bufs[10] + ITER_PREV, bufs[11] + ITER_PREV,
				bufs[12] + ITER_PREV, bufs[13] + ITER_PREV, bufs[14] + ITER_PREV, bufs[15] + ITER_PREV,
				h,
				bufs[0], bufs[1], bufs[2], bufs[3], bufs[4], bufs[5], bufs[6], bufs[7],
				bufs[8], bufs[9], bufs[10], bufs[11], bufs[12], bufs[13], bufs[14], bufs[15],
				ITER_LEN(n));
		}
	}
	for (size_t r = 0; r < rem; r++)
		bufs[r] = buffer[local + offset - rem + r];
	for (size_t r = rem; r < m_num; r++)
		bufs[r] = &empty[0];
	for (size_t j = 0; j < (1 << w) - 1; j++) {
		for (size_t r = 0; r < rem; r++)
			bufs[r][ITER_J] = j;
		hss_hash_ctx_sixteen(
			bufs[0] + ITER_PREV, bufs[1] + ITER_PREV, bufs[2] + ITER_PREV, bufs[3] + ITER_PREV,
			bufs[4] + ITER_PREV, bufs[5] + ITER_PREV, bufs[6] + ITER_PREV, bufs[7] + ITER_PREV,
			bufs[8] + ITER_PREV, bufs[9] + ITER_PREV, bufs[10] + ITER_PREV, bufs[11] + ITER_PREV,
			bufs[12] + ITER_PREV, bufs[13] + ITER_PREV, bufs[14] + ITER_PREV, bufs[15] + ITER_PREV,
			h,
			bufs[0], bufs[1], bufs[2], bufs[3], bufs[4], bufs[5], bufs[6], bufs[7],
			bufs[8], bufs[9], bufs[10], bufs[11], bufs[12], bufs[13], bufs[14], bufs[15],
			ITER_LEN(n));
	}
	//通信

	// end = MPI_Wtime();
	// if(me == 0) printf("1 %lf\n", 1000 * (end - start));
	// start = MPI_Wtime();
#ifdef Gather_communicator
	int r_counts[new_size];
	int r_dis[new_size];        //存起始位置

	for (int i = 0; i < f; i++) {
		r_counts[i] = (e + 1) * ITER_MAX_LEN;
		r_dis[i] = (i * e + i) * ITER_MAX_LEN;
	}
	for (int i = f; i < new_size; i++) {
		r_counts[i] = e * ITER_MAX_LEN;
		r_dis[i] = (i * e + f) * ITER_MAX_LEN;
	}
	MPI_Gatherv(buffer[offset], local * ITER_MAX_LEN, MPI_CHAR, buffer, r_counts, r_dis, MPI_CHAR, 0, local_comm);
#else /*思考：是否没有划通信子域，不可以用MPI_Gatherv？，选择使用send和recv*/
	// printf("************00me = %d\n", me);
	int a = me % global_level_nodes;
	if (id != 0)
		MPI_Send(buffer[offset], local * ITER_MAX_LEN, MPI_CHAR, a, 0, MPI_COMM_WORLD);
	if (id == 0) {
		if (f > 1) {
			for (size_t i = 1; i < f; i++) {
				unsigned int b = a + (i * global_level_nodes);
				MPI_Recv(buffer[i * e + i], (e + 1) * ITER_MAX_LEN, MPI_CHAR, b, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			}
		}
		for (size_t i = f; i < new_size; i++) {
			unsigned int b = a + (i * global_level_nodes);
			MPI_Recv(buffer[i * e + f], e * ITER_MAX_LEN, MPI_CHAR, b, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
	}

#endif

	for (size_t i = 0; i < p; i++)
		hss_update_hash_context(h, &public_ctx, buffer[i] + ITER_PREV, n);
	/* And the result of the running hash is the public key */
	hss_finalize_hash_context(h, &public_ctx, public_key);

	for (size_t i = 0; i < 16; i++)
		hss_zeroize(&ctx[i], sizeof ctx[i]);
	// end = MPI_Wtime();
	// if(me == 0) printf("2 %lf\n", 1000 * (end - start));

	return true;
}
bool mpi_sign_lm_ots_generate_public_key(
	param_set_t lm_ots_type,
	const unsigned char *I, /* Public key identifier */
	merkle_index_t q,       /* Diversification string, 4 bytes value */
	struct seed_derive *seed,
	unsigned char *public_key, size_t public_key_len)
{
	int me, proc_num;

	MPI_Comm_rank(MPI_COMM_WORLD, &me);             //进程号
	MPI_Comm_size(MPI_COMM_WORLD, &proc_num);       //进程数
	// double start, end;

	// start = MPI_Wtime();

	// end = MPI_Wtime();
	// if(me == 0) printf("0 %lf\n", 1000 * (end - start));
	// start = MPI_Wtime();

	/* Look up the parameter set */
	unsigned h, n, w, p, ls;

	if (!lm_ots_look_up_parameter_set(lm_ots_type, &h, &n, &w, &p, &ls))
		return false;

	/* Start the hash that computes the final value */
	union hash_context public_ctx;

	hss_init_hash_context(h, &public_ctx);
	{
		unsigned char prehash_prefix[PBLC_PREFIX_LEN];
		memcpy(prehash_prefix + PBLC_I, I, I_LEN);
		put_bigendian(prehash_prefix + PBLC_Q, q, 4);
		SET_D(prehash_prefix + PBLC_D, D_PBLC);
		hss_update_hash_context(h, &public_ctx, prehash_prefix,
					PBLC_PREFIX_LEN);
	}
	int id, new_size;                  //子通信域中新的进程号

#ifdef Gather_communicator

	MPI_Comm_rank(local_comm02, &id);               //进程号
	MPI_Comm_size(local_comm02, &new_size);         //进程数
#else
	id = me / num01;
	new_size = proc_num / num01;
#endif
	// int i, j;
	// new_size是子通信域的总进程数
	int e = p / new_size;
	int f = p % new_size;
	unsigned int local = e + ((id < f) ? 1 : 0);
	unsigned int offset = id * e + ((id < f) ? id : f);

	hss_seed_derive_set_j(seed, 0);
	int num = local / 16;
	int rem = local % 16;
	int m_num = 16;
	unsigned char buffer[((p / 16) + 1) * 16][ITER_MAX_LEN];
	unsigned char *bufs[m_num];
	union hash_context ctx[16];
	unsigned char empty[ITER_MAX_LEN];

	for (size_t i = 0; i < m_num; i++)
		bufs[i] = &empty[0];

	for (size_t i = 0; i < p; i++) {
		memcpy(buffer[i] + ITER_I, I, I_LEN);
		put_bigendian(buffer[i] + ITER_Q, q, 4);
		hss_seed_derive(buffer[i] + ITER_PREV, seed, i < p - 1);
		put_bigendian(buffer[i] + ITER_K, i, 2);
	}
	for (size_t i = 0; i < num; i++) {
		int suboffset = offset + i * 16;
		for (size_t r = 0; r < m_num; r++)
			bufs[r] = buffer[suboffset + r];
		for (size_t j = 0; j < (1 << w) - 1; j++) {
			bufs[0][ITER_J] = j;
			bufs[1][ITER_J] = j;
			bufs[2][ITER_J] = j;
			bufs[3][ITER_J] = j;
			bufs[4][ITER_J] = j;
			bufs[5][ITER_J] = j;
			bufs[6][ITER_J] = j;
			bufs[7][ITER_J] = j;
			bufs[8][ITER_J] = j;
			bufs[9][ITER_J] = j;
			bufs[10][ITER_J] = j;
			bufs[11][ITER_J] = j;
			bufs[12][ITER_J] = j;
			bufs[13][ITER_J] = j;
			bufs[14][ITER_J] = j;
			bufs[15][ITER_J] = j;
			hss_hash_ctx_sixteen(
				bufs[0] + ITER_PREV, bufs[1] + ITER_PREV, bufs[2] + ITER_PREV, bufs[3] + ITER_PREV,
				bufs[4] + ITER_PREV, bufs[5] + ITER_PREV, bufs[6] + ITER_PREV, bufs[7] + ITER_PREV,
				bufs[8] + ITER_PREV, bufs[9] + ITER_PREV, bufs[10] + ITER_PREV, bufs[11] + ITER_PREV,
				bufs[12] + ITER_PREV, bufs[13] + ITER_PREV, bufs[14] + ITER_PREV, bufs[15] + ITER_PREV,
				h,
				bufs[0], bufs[1], bufs[2], bufs[3], bufs[4], bufs[5], bufs[6], bufs[7],
				bufs[8], bufs[9], bufs[10], bufs[11], bufs[12], bufs[13], bufs[14], bufs[15],
				ITER_LEN(n));
		}
	}
	for (size_t r = 0; r < rem; r++)
		bufs[r] = buffer[local + offset - rem + r];
	for (size_t r = rem; r < m_num; r++)
		bufs[r] = &empty[0];
	for (size_t j = 0; j < (1 << w) - 1; j++) {
		for (size_t r = 0; r < rem; r++)
			bufs[r][ITER_J] = j;
		hss_hash_ctx_sixteen(
			bufs[0] + ITER_PREV, bufs[1] + ITER_PREV, bufs[2] + ITER_PREV, bufs[3] + ITER_PREV,
			bufs[4] + ITER_PREV, bufs[5] + ITER_PREV, bufs[6] + ITER_PREV, bufs[7] + ITER_PREV,
			bufs[8] + ITER_PREV, bufs[9] + ITER_PREV, bufs[10] + ITER_PREV, bufs[11] + ITER_PREV,
			bufs[12] + ITER_PREV, bufs[13] + ITER_PREV, bufs[14] + ITER_PREV, bufs[15] + ITER_PREV,
			h,
			bufs[0], bufs[1], bufs[2], bufs[3], bufs[4], bufs[5], bufs[6], bufs[7],
			bufs[8], bufs[9], bufs[10], bufs[11], bufs[12], bufs[13], bufs[14], bufs[15],
			ITER_LEN(n));
	}
#ifdef Gather_communicator
	//通信
	int r_counts[new_size];
	int r_dis[new_size];        //存起始位置

	for (int i = 0; i < f; i++) {
		r_counts[i] = (e + 1) * ITER_MAX_LEN;
		r_dis[i] = (i * e + i) * ITER_MAX_LEN;
	}
	for (int i = f; i < new_size; i++) {
		r_counts[i] = e * ITER_MAX_LEN;
		r_dis[i] = (i * e + f) * ITER_MAX_LEN;
	}
	// end = MPI_Wtime();
	// if(me == 0) printf("1 %lf\n", 1000 * (end - start));
	// start = MPI_Wtime();

	MPI_Gatherv(buffer[offset], local * ITER_MAX_LEN, MPI_CHAR, buffer, r_counts, r_dis, MPI_CHAR, 0, local_comm02);
	if (id == 0) {
		for (size_t i = 0; i < p; i++)
			hss_update_hash_context(h, &public_ctx, buffer[i] + ITER_PREV, n);
		/* And the result of the running hash is the public key */
		hss_finalize_hash_context(h, &public_ctx, public_key);
	}
	MPI_Bcast(public_key, public_key_len, MPI_CHAR, 0, local_comm02);


#else

	int a = me % num01;
	if (id != 0)
		// printf("**************02 me = %d\n", me);
		MPI_Send(buffer[offset], local * ITER_MAX_LEN, MPI_CHAR, a, 0, MPI_COMM_WORLD);
	if (id == 0) {
		if (f > 1) {
			for (size_t i = 1; i < f; i++) {
				unsigned int b = a + (i * num01);
				MPI_Recv(buffer[i * e + i], (e + 1) * ITER_MAX_LEN, MPI_CHAR, b, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				// printf("**************03 me = %d\n", me);
			}
		}
		for (size_t i = f; i < new_size; i++) {
			unsigned int b = a + (i * num01);
			MPI_Recv(buffer[i * e + f], e * ITER_MAX_LEN, MPI_CHAR, b, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			// printf("**************04 me = %d\n", me);
		}
	}
	if (id == 0) {
		for (size_t i = 0; i < p; i++)
			hss_update_hash_context(h, &public_ctx, buffer[i] + ITER_PREV, n);
		/* And the result of the running hash is the public key */
		hss_finalize_hash_context(h, &public_ctx, public_key);


		for (size_t i = 1; i < new_size; i++) {
			unsigned a = me + i * num01;
			MPI_Send(public_key, public_key_len, MPI_CHAR, a, 0, MPI_COMM_WORLD);
		}
	}
	if (id != 0) {
		unsigned b = (me % num01);
		MPI_Recv(public_key, public_key_len, MPI_CHAR, b, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}

#endif

	// MPI_Bcast(public_key, public_key_len, MPI_CHAR, 0, local_comm02);
	// printf("**************03 me = %d\n", me);
	for (size_t i = 0; i < 16; i++)
		hss_zeroize(&ctx[i], sizeof ctx[i]);

	// end = MPI_Wtime();
	// if(me == 0) printf("2 %lf\n", 1000 * (end - start));

	return true;
}

bool lm_ots_generate_public_key(
	param_set_t lm_ots_type,
	const unsigned char *I, /* Public key identifier */
	merkle_index_t q,       /* Diversification string, 4 bytes value */
	struct seed_derive *seed,
	unsigned char *public_key, size_t public_key_len)
{
	/* Look up the parameter set */
	unsigned h, n, w, p, ls;
	if (!lm_ots_look_up_parameter_set(lm_ots_type, &h, &n, &w, &p, &ls))
		return false;

	/* Start the hash that computes the final value */
	union hash_context public_ctx;

	hss_init_hash_context(h, &public_ctx);
	{
		unsigned char prehash_prefix[PBLC_PREFIX_LEN];
		memcpy(prehash_prefix + PBLC_I, I, I_LEN);
		put_bigendian(prehash_prefix + PBLC_Q, q, 4);
		SET_D(prehash_prefix + PBLC_D, D_PBLC);
		hss_update_hash_context(h, &public_ctx, prehash_prefix,
					PBLC_PREFIX_LEN);
	}

	int i, j;

	hss_seed_derive_set_j(seed, 0);
	int num = p / 16;
	int rem = p % 16;
	unsigned char buffer[(num + 1) * 16][ITER_MAX_LEN];
	union hash_context ctx[16];

	for (size_t i = 0; i < p; i++) {
		memcpy(buffer[i] + ITER_I, I, I_LEN);
		put_bigendian(buffer[i] + ITER_Q, q, 4);
		hss_seed_derive(buffer[i] + ITER_PREV, seed, i < p - 1);
		put_bigendian(buffer[i] + ITER_K, i, 2);
	}
	for (size_t i = 0; i < num; i++) {
		int offset = i * 16;
		for (j = 0; j < (1 << w) - 1; j++) {
			buffer[offset][ITER_J] = j;
			buffer[offset + 1][ITER_J] = j;
			buffer[offset + 2][ITER_J] = j;
			buffer[offset + 3][ITER_J] = j;
			buffer[offset + 4][ITER_J] = j;
			buffer[offset + 5][ITER_J] = j;
			buffer[offset + 6][ITER_J] = j;
			buffer[offset + 7][ITER_J] = j;
			buffer[offset + 8][ITER_J] = j;
			buffer[offset + 9][ITER_J] = j;
			buffer[offset + 10][ITER_J] = j;
			buffer[offset + 11][ITER_J] = j;
			buffer[offset + 12][ITER_J] = j;
			buffer[offset + 13][ITER_J] = j;
			buffer[offset + 14][ITER_J] = j;
			buffer[offset + 15][ITER_J] = j;

			hss_hash_ctx_sixteen(buffer[offset] + ITER_PREV, buffer[offset + 1] + ITER_PREV, buffer[offset + 2] + ITER_PREV, buffer[offset + 3] + ITER_PREV, buffer[offset + 4] + ITER_PREV, buffer[offset + 5] + ITER_PREV, buffer[offset + 6] + ITER_PREV, buffer[offset + 7] + ITER_PREV, buffer[offset + 8] + ITER_PREV, buffer[offset + 9] + ITER_PREV, buffer[offset + 10] + ITER_PREV, buffer[offset + 11] + ITER_PREV, buffer[offset + 12] + ITER_PREV, buffer[offset + 13] + ITER_PREV, buffer[offset + 14] + ITER_PREV, buffer[offset + 15] + ITER_PREV, h, buffer[offset], buffer[offset + 1], buffer[offset + 2], buffer[offset + 3], buffer[offset + 4], buffer[offset + 5], buffer[offset + 6], buffer[offset + 7], buffer[offset + 8], buffer[offset + 9], buffer[offset + 10], buffer[offset + 11], buffer[offset + 12], buffer[offset + 13], buffer[offset + 14], buffer[offset + 15], ITER_LEN(n));
		}
	}
	int m = p - rem;

	for (j = 0; j < (1 << w) - 1; j++) {
		for (size_t m = p - rem; m < p; m++)
			buffer[m][ITER_J] = j;
		m = p - rem;

		hss_hash_ctx_sixteen(buffer[m] + ITER_PREV, buffer[m + 1] + ITER_PREV, buffer[m + 2] + ITER_PREV, buffer[m + 3] + ITER_PREV, buffer[m + 4] + ITER_PREV, buffer[m + 5] + ITER_PREV, buffer[m + 6] + ITER_PREV, buffer[m + 7] + ITER_PREV, buffer[m + 8] + ITER_PREV, buffer[m + 9] + ITER_PREV, buffer[m + 10] + ITER_PREV, buffer[m + 11] + ITER_PREV, buffer[m + 12] + ITER_PREV, buffer[m + 13] + ITER_PREV, buffer[m + 14] + ITER_PREV, buffer[m + 15] + ITER_PREV, h, buffer[m], buffer[m + 1], buffer[m + 2], buffer[m + 3], buffer[m + 4], buffer[m + 5], buffer[m + 6], buffer[m + 7], buffer[m + 8], buffer[m + 9], buffer[m + 10], buffer[m + 11], buffer[m + 12], buffer[m + 13], buffer[m + 14], buffer[m + 15], ITER_LEN(n));
	}

	for (size_t i = 0; i < p; i++)
		hss_update_hash_context(h, &public_ctx, buffer[i] + ITER_PREV, n);
	/* And the result of the running hash is the public key */
	hss_finalize_hash_context(h, &public_ctx, public_key);

	for (size_t i = 0; i < 16; i++)
		hss_zeroize(&ctx[i], sizeof ctx[i]);
	return true;
}

#else
/*vectorization -- eight */
/*keygen version*/
bool mpi_lm_ots_generate_public_key(
	param_set_t lm_ots_type,
	const unsigned char *I, /* Public key identifier */
	merkle_index_t q,       /* Diversification string, 4 bytes value */
	struct seed_derive *seed,
	unsigned char *public_key, size_t public_key_len)
{
	int me, proc_num;

	MPI_Comm_rank(MPI_COMM_WORLD, &me);             //进程号
	MPI_Comm_size(MPI_COMM_WORLD, &proc_num);       //进程数
	// double start, end;

	// start = MPI_Wtime();
	unsigned int new_size;

	new_size = proc_num / global_level_nodes;

	// new_size = 2;
	MPI_Barrier(MPI_COMM_WORLD);

	// end = MPI_Wtime();
	// if(me == 0) printf("0 %lf\n", 1000 * (end - start));
	// start = MPI_Wtime();

	/* Look up the parameter set */
	unsigned h, n, w, p, ls;

	if (!lm_ots_look_up_parameter_set(lm_ots_type, &h, &n, &w, &p, &ls))
		return false;
	MPI_Barrier(MPI_COMM_WORLD);
	/* Start the hash that computes the final value */
	union hash_context public_ctx;

	hss_init_hash_context(h, &public_ctx);
	{
		unsigned char prehash_prefix[PBLC_PREFIX_LEN];
		memcpy(prehash_prefix + PBLC_I, I, I_LEN);
		put_bigendian(prehash_prefix + PBLC_Q, q, 4);
		SET_D(prehash_prefix + PBLC_D, D_PBLC);
		hss_update_hash_context(h, &public_ctx, prehash_prefix,
					PBLC_PREFIX_LEN);
	}
	int id;                                 //子通信域中新的进程号

#if (defined(Gather02)) || (defined(Igather02))
	MPI_Comm_rank(local_comm, &id);         //进程号
#else
	//假如没有划分通信子域
	id = me / global_level_nodes;
#endif
	// new_size是子通信域的总进程数
	int e = p / new_size;
	int f = p % new_size;
	unsigned int local = e + ((id < f) ? 1 : 0);
	unsigned int offset = id * e + ((id < f) ? id : f);

	hss_seed_derive_set_j(seed, 0);
	int num = local / 8;
	int rem = local % 8;
	int m_num = 8;
	unsigned char buffer[((p / 8) + 1) * 8][ITER_MAX_LEN];
	unsigned char *bufs[m_num];
	union hash_context ctx[8];
	unsigned char empty[ITER_MAX_LEN];

	for (size_t i = 0; i < 72; i++)
		bufs[i] = &empty[0];

	for (size_t i = 0; i < p; i++) {
		memcpy(buffer[i] + ITER_I, I, I_LEN);
		put_bigendian(buffer[i] + ITER_Q, q, 4);
		hss_seed_derive(buffer[i] + ITER_PREV, seed, i < p - 1);
		put_bigendian(buffer[i] + ITER_K, i, 2);
	}
	for (size_t i = 0; i < num; i++) {
		int suboffset = offset + i * 8;
		for (size_t r = 0; r < m_num; r++)
			bufs[r] = buffer[suboffset + r];
		for (size_t j = 0; j < (1 << w) - 1; j++) {
			bufs[0][ITER_J] = j;
			bufs[1][ITER_J] = j;
			bufs[2][ITER_J] = j;
			bufs[3][ITER_J] = j;
			bufs[4][ITER_J] = j;
			bufs[5][ITER_J] = j;
			bufs[6][ITER_J] = j;
			bufs[7][ITER_J] = j;
			hss_hash_ctx_eight(
				bufs[0] + ITER_PREV, bufs[1] + ITER_PREV, bufs[2] + ITER_PREV, bufs[3] + ITER_PREV,
				bufs[4] + ITER_PREV, bufs[5] + ITER_PREV, bufs[6] + ITER_PREV, bufs[7] + ITER_PREV,
				h,
				bufs[0], bufs[1], bufs[2], bufs[3], bufs[4], bufs[5], bufs[6], bufs[7], ITER_LEN(n));
		}
	}
	for (size_t r = 0; r < rem; r++)
		bufs[r] = buffer[local + offset - rem + r];
	for (size_t r = rem; r < m_num; r++)
		bufs[r] = &empty[0];
	for (size_t j = 0; j < (1 << w) - 1; j++) {
		for (size_t r = 0; r < rem; r++)
			bufs[r][ITER_J] = j;
		hss_hash_ctx_eight(
			bufs[0] + ITER_PREV, bufs[1] + ITER_PREV, bufs[2] + ITER_PREV, bufs[3] + ITER_PREV,
			bufs[4] + ITER_PREV, bufs[5] + ITER_PREV, bufs[6] + ITER_PREV, bufs[7] + ITER_PREV,
			h,
			bufs[0], bufs[1], bufs[2], bufs[3], bufs[4], bufs[5], bufs[6], bufs[7], ITER_LEN(n));
	}
	//通信

	// end = MPI_Wtime();
	// if(me == 0) printf("1 %lf\n", 1000 * (end - start));
	// start = MPI_Wtime();

#if (defined(Gather02))
	int r_counts[new_size];
	int r_dis[new_size];        //存起始位置

	for (int i = 0; i < f; i++) {
		r_counts[i] = (e + 1) * ITER_MAX_LEN;
		r_dis[i] = (i * e + i) * ITER_MAX_LEN;
	}
	for (int i = f; i < new_size; i++) {
		r_counts[i] = e * ITER_MAX_LEN;
		r_dis[i] = (i * e + f) * ITER_MAX_LEN;
	}
#if (defined(Communication_time))
	MPI_Barrier(MPI_COMM_WORLD); // 同步
	double t0, t1, t2, t_sum;
	t0 = MPI_Wtime();
#endif
	MPI_Gatherv(buffer[offset], local * ITER_MAX_LEN, MPI_CHAR, buffer, r_counts, r_dis, MPI_CHAR, 0, local_comm);

#if (defined(Communication_time))
	t1 = MPI_Wtime();
	t2 = t1 - t0;
	// printf("t2 = %f\n", t2 * 1000);

	MPI_Reduce(&t2, &t_sum, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	sum_time2 += (t_sum * 1000);
	if (me == 0)
		printf("sum LM-OTS_public_key Gather_time = %.6lf\n", sum_time2);
#endif
#elif (defined(Igather02))
	int r_counts[new_size];
	int r_dis[new_size]; //存起始位置

	for (int i = 0; i < f; i++) {
		r_counts[i] = (e + 1) * ITER_MAX_LEN;
		r_dis[i] = (i * e + i) * ITER_MAX_LEN;
	}
	for (int i = f; i < new_size; i++) {
		r_counts[i] = e * ITER_MAX_LEN;
		r_dis[i] = (i * e + f) * ITER_MAX_LEN;
	}
#if (defined(Communication_time))
	MPI_Barrier(MPI_COMM_WORLD); // 同步
	double t0, t1, t2, t_sum;
	t0 = MPI_Wtime();
#endif
	MPI_Request re;
	MPI_Igatherv(buffer[offset], local * ITER_MAX_LEN, MPI_CHAR, buffer, r_counts, r_dis, MPI_CHAR, 0, local_comm, &re);
	MPI_Wait(&re, MPI_STATUS_IGNORE);
#if (defined(Communication_time))
	t1 = MPI_Wtime();
	t2 = t1 - t0;
	// printf("t2 = %f\n", t2 * 1000);

	MPI_Reduce(&t2, &t_sum, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	sum_time2 += (t_sum * 1000);
	if (me == 0)
		printf("sum LM-OTS_public_key Igather_time = %.6lf\n", sum_time2);
#endif
#elif (defined(Send_Recv02))
	int a = me % global_level_nodes;
	MPI_Barrier(MPI_COMM_WORLD); // 同步
	double t0, t1, t2, t_sum;
	t0 = MPI_Wtime();
	if (id != 0)
		MPI_Send(buffer[offset], local * ITER_MAX_LEN, MPI_CHAR, a, 0, MPI_COMM_WORLD);
	if (id == 0) {
		if (f > 1) {
			for (size_t i = 1; i < f; i++) {
				unsigned int b = a + (i * global_level_nodes);
				MPI_Recv(buffer[i * e + i], (e + 1) * ITER_MAX_LEN, MPI_CHAR, b, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			}
		}
		for (size_t i = f; i < new_size; i++) {
			unsigned int b = a + (i * global_level_nodes);
			MPI_Recv(buffer[i * e + f], e * ITER_MAX_LEN, MPI_CHAR, b, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
	}
#if (defined(Communication_time))
	t1 = MPI_Wtime();
	t2 = t1 - t0;
	// printf("t2 = %f\n", t2 * 1000);

	MPI_Reduce(&t2, &t_sum, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	sum_time2 += (t_sum * 1000);
	if (me == 0)
		printf("sum LM-OTS_public_key Send_Recv_time = %.6lf\n", sum_time2);
#endif
#elif (defined(Isend_Irecv02))
/*非阻塞通信*/
	int a = me % global_level_nodes;
	MPI_Request request1, request2[new_size], request3[new_size], request4;
#if (defined(Communication_time))
	MPI_Barrier(MPI_COMM_WORLD); // 同步
	double t0, t1, t2, t_sum;
	t0 = MPI_Wtime();
#endif
	if (id != 0)
		MPI_Isend(buffer[offset], local * ITER_MAX_LEN, MPI_CHAR, a, 0, MPI_COMM_WORLD, &request1);
	if (id == 0) {
		if (f > 1) {
			for (size_t i = 1; i < f; i++) {
				unsigned int b = a + (i * global_level_nodes);
				MPI_Irecv(buffer[i * e + i], (e + 1) * ITER_MAX_LEN, MPI_CHAR, b, 0, MPI_COMM_WORLD, &request2[i]);
			}
		}
		for (size_t i = f; i < new_size; i++) {
			unsigned int b = a + (i * global_level_nodes);
			MPI_Irecv(buffer[i * e + f], e * ITER_MAX_LEN, MPI_CHAR, b, 0, MPI_COMM_WORLD, &request2[i]);
		}
	}
	if (id != 0)
		MPI_Waitall(1, &request1, MPI_STATUS_IGNORE);
	if (id == 0) {
		if (f > 1)
			MPI_Waitall(f - 1, &request2[1], MPI_STATUS_IGNORE);
		MPI_Waitall(new_size - f, &request2[f], MPI_STATUS_IGNORE);
	}

	// if (id == 0) {
	// 	for (size_t i = 0; i < new_size; i++) {
	// 		unsigned int b = me + (i * global_level_nodes);
	// 		MPI_Isend(buffer, p * ITER_MAX_LEN, MPI_CHAR, b, 0, MPI_COMM_WORLD, &request3[i]);
	// 	}
	// }
	// if (id != 0)
	// 	// unsigned int b = a + (i * global_level_nodes);
	// 	MPI_Irecv(buffer, p * ITER_MAX_LEN, MPI_CHAR, a, 0, MPI_COMM_WORLD, &request4);
	// if (id == 0) {
	// 	MPI_Waitall(new_size, &request3, MPI_STATUS_IGNORE);
	// }
	// if (id != 0)
	// 	MPI_Waitall(1, &request4, MPI_STATUS_IGNORE);
#if (defined(Communication_time))
	t1 = MPI_Wtime();
	t2 = t1 - t0;
	// printf("t2 = %f\n", t2 * 1000);

	MPI_Reduce(&t2, &t_sum, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	sum_time2 += (t_sum * 1000);
	if (me == 0)
		printf("************MAX LM-OTS_public_key Isend_Irecv Communication_time = %.6lf***************\n", sum_time2);
#endif
#endif

	for (size_t i = 0; i < p; i++)
		hss_update_hash_context(h, &public_ctx, buffer[i] + ITER_PREV, n);
	/* And the result of the running hash is the public key */
	hss_finalize_hash_context(h, &public_ctx, public_key);

	for (size_t i = 0; i < 8; i++)
		hss_zeroize(&ctx[i], sizeof ctx[i]);
	// end = MPI_Wtime();
	// if(me == 0) printf("2 %lf\n", 1000 * (end - start));

	return true;
}

bool mpi_sign_lm_ots_generate_public_key(
	param_set_t lm_ots_type,
	const unsigned char *I, /* Public key identifier */
	merkle_index_t q,       /* Diversification string, 4 bytes value */
	struct seed_derive *seed,
	unsigned char *public_key, size_t public_key_len)
{
	int me, proc_num;

	MPI_Comm_rank(MPI_COMM_WORLD, &me);             //进程号
	MPI_Comm_size(MPI_COMM_WORLD, &proc_num);       //进程数
	// double start, end;

	// start = MPI_Wtime();

	// end = MPI_Wtime();
	// if(me == 0) printf("0 %lf\n", 1000 * (end - start));
	// start = MPI_Wtime();

	/* Look up the parameter set */
	unsigned h, n, w, p, ls;

	if (!lm_ots_look_up_parameter_set(lm_ots_type, &h, &n, &w, &p, &ls))
		return false;

	/* Start the hash that computes the final value */
	union hash_context public_ctx;

	hss_init_hash_context(h, &public_ctx);
	{
		unsigned char prehash_prefix[PBLC_PREFIX_LEN];
		memcpy(prehash_prefix + PBLC_I, I, I_LEN);
		put_bigendian(prehash_prefix + PBLC_Q, q, 4);
		SET_D(prehash_prefix + PBLC_D, D_PBLC);
		hss_update_hash_context(h, &public_ctx, prehash_prefix,
					PBLC_PREFIX_LEN);
	}
	int id, new_size;                  //子通信域中新的进程号

#ifdef Gather06

	MPI_Comm_rank(local_comm02, &id);               //进程号
	MPI_Comm_size(local_comm02, &new_size);         //进程数
#else
	id = me / num01;
	new_size = proc_num / num01;
#endif
	int max_pro = (p / 8) + 1;
	int work_pro = 0;
	if (new_size > max_pro) work_pro = max_pro;
	else work_pro = new_size;
	// if(id == 0) printf("work_pro %d\n", work_pro);
	int e = p / work_pro;
	int f = p % work_pro;
	unsigned int local = e + ((id < f) ? 1 : 0);
	unsigned int offset = id * e + ((id < f) ? id : f);

	hss_seed_derive_set_j(seed, 0);
	int num = local / 8;
	int rem = local % 8;
	int m_num = 8;
	unsigned char buffer[((p / 8) + 1) * 8][ITER_MAX_LEN];
	unsigned char *bufs[m_num];
	union hash_context ctx[8];
	unsigned char empty[ITER_MAX_LEN];
	if (id < work_pro){

		for (size_t i = 0; i < 72; i++)
		bufs[i] = &empty[0];

		for (size_t i = 0; i < p; i++) {
			memcpy(buffer[i] + ITER_I, I, I_LEN);
			put_bigendian(buffer[i] + ITER_Q, q, 4);
			hss_seed_derive(buffer[i] + ITER_PREV, seed, i < p - 1);
			put_bigendian(buffer[i] + ITER_K, i, 2);
		}
		for (size_t i = 0; i < num; i++) {
			int suboffset = offset + i * 8;
			for (size_t r = 0; r < m_num; r++)
			bufs[r] = buffer[suboffset + r];
			for (size_t j = 0; j < (1 << w) - 1; j++) {
				bufs[0][ITER_J] = j;
				bufs[1][ITER_J] = j;
				bufs[2][ITER_J] = j;
				bufs[3][ITER_J] = j;
				bufs[4][ITER_J] = j;
				bufs[5][ITER_J] = j;
				bufs[6][ITER_J] = j;
				bufs[7][ITER_J] = j;
				hss_hash_ctx_eight(
					bufs[0] + ITER_PREV, bufs[1] + ITER_PREV, bufs[2] + ITER_PREV, bufs[3] + ITER_PREV,
					bufs[4] + ITER_PREV, bufs[5] + ITER_PREV, bufs[6] + ITER_PREV, bufs[7] + ITER_PREV,
					h,
					bufs[0], bufs[1], bufs[2], bufs[3], bufs[4], bufs[5], bufs[6], bufs[7], ITER_LEN(n));
				}
			}
			for (size_t r = 0; r < rem; r++)
			bufs[r] = buffer[local + offset - rem + r];
			for (size_t r = rem; r < m_num; r++)
			bufs[r] = &empty[0];
			for (size_t j = 0; j < (1 << w) - 1; j++) {
				for (size_t r = 0; r < rem; r++)
				bufs[r][ITER_J] = j;
				hss_hash_ctx_eight(
					bufs[0] + ITER_PREV, bufs[1] + ITER_PREV, bufs[2] + ITER_PREV, bufs[3] + ITER_PREV,
					bufs[4] + ITER_PREV, bufs[5] + ITER_PREV, bufs[6] + ITER_PREV, bufs[7] + ITER_PREV,
					h,
					bufs[0], bufs[1], bufs[2], bufs[3], bufs[4], bufs[5], bufs[6], bufs[7], ITER_LEN(n));
				}
	}
#if (defined(Gather06))
	//通信
	int r_counts[new_size];
	int r_dis[new_size];        //存起始位置

	for (int i = 0; i < f; i++) {
		r_counts[i] = (e + 1) * ITER_MAX_LEN;
		r_dis[i] = (i * e + i) * ITER_MAX_LEN;
	}
	for (int i = f; i < new_size; i++) {
		r_counts[i] = e * ITER_MAX_LEN;
		r_dis[i] = (i * e + f) * ITER_MAX_LEN;
	}
	// end = MPI_Wtime();
	// if(me == 0) printf("1 %lf\n", 1000 * (end - start));
	// start = MPI_Wtime();

	//划分通信子域
	MPI_Comm local_comm_Gather06;
	if (new_size > work_pro){
	    int sub_num; 	//子通信域中的进程数
	    sub_num = work_pro;
	    int local_rank_Gather06[work_pro];
	    MPI_Group World_Group, local_group;
	    MPI_Comm_group(local_comm02, &World_Group);
	    //8是level_nodes
	    for (int i = 0; i < work_pro; i++) {
	        local_rank_Gather06[i] = i;         //存local_comm02子域进程号
	    }
	    MPI_Group_incl(World_Group, sub_num, local_rank_Gather06, &local_group);
	    MPI_Comm_create(MPI_COMM_WORLD, local_group, &local_comm_Gather06);

	}

	if (new_size > work_pro){
	    if (id < work_pro){
	        MPI_Gatherv(buffer[offset], local * ITER_MAX_LEN, MPI_CHAR, buffer, r_counts, r_dis, MPI_CHAR, 0, local_comm_Gather06);
	    }
	}
	else
		MPI_Gatherv(buffer[offset], local * ITER_MAX_LEN, MPI_CHAR, buffer, r_counts, r_dis, MPI_CHAR, 0, local_comm02);
	if (id == 0) {
		for (size_t i = 0; i < p; i++)
			hss_update_hash_context(h, &public_ctx, buffer[i] + ITER_PREV, n);
		/* And the result of the running hash is the public key */
		hss_finalize_hash_context(h, &public_ctx, public_key);
	}
	MPI_Bcast(public_key, public_key_len, MPI_CHAR, 0, local_comm02);


#elif (defined(Send_Recv06))

	int a = me % num01;
	if (id != 0  && id < work_pro)
		// printf("**************02 me = %d\n", me);
		MPI_Send(buffer[offset], local * ITER_MAX_LEN, MPI_CHAR, a, 0, MPI_COMM_WORLD);
	if (id == 0) {
		if (f > 1) {
			for (size_t i = 1; i < f; i++) {
				unsigned int b = a + (i * num01);
				MPI_Recv(buffer[i * e + i], (e + 1) * ITER_MAX_LEN, MPI_CHAR, b, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				// printf("**************03 me = %d\n", me);
			}
			for (size_t i = f; i < work_pro; i++) {
				unsigned int b = a + (i * num01);
				MPI_Recv(buffer[i * e + f], e * ITER_MAX_LEN, MPI_CHAR, b, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				// printf("**************04 me = %d\n", me);
			}
		} else {
			for (size_t i = 1; i < work_pro; i++) {
				unsigned int b = a + (i * num01);
				MPI_Recv(buffer[i * e + f], e * ITER_MAX_LEN, MPI_CHAR, b, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				// printf("**************04 me = %d\n", me);
			}
		}
	}

	// if(me == 0)
	// printf("me = %d, %d, %d\n", num01, me, f);

	if (id == 0) {
		for (size_t i = 0; i < p; i++)
			hss_update_hash_context(h, &public_ctx, buffer[i] + ITER_PREV, n);
		/* And the result of the running hash is the public key */
		hss_finalize_hash_context(h, &public_ctx, public_key);


		for (size_t i = 1; i < new_size; i++) {
			unsigned a = me + i * num01;
			MPI_Send(public_key, public_key_len, MPI_CHAR, a, 0, MPI_COMM_WORLD);
		}
	}
	if (id != 0) {
		unsigned b = (me % num01);
		MPI_Recv(public_key, public_key_len, MPI_CHAR, b, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}

#endif
	// if (id == 0) {
	// 	for (size_t i = 0; i < p; i++)
	// 		hss_update_hash_context(h, &public_ctx, buffer[i] + ITER_PREV, n);
	// 	/* And the result of the running hash is the public key */
	// 	hss_finalize_hash_context(h, &public_ctx, public_key);
	//
	//
	// 	for (size_t i = 1; i < new_size; i++) {
	// 		unsigned a = me + i * num01;
	// 		MPI_Send(public_key, public_key_len, MPI_CHAR, a, 0, MPI_COMM_WORLD);
	// 	}
	// }
	// if (id != 0) {
	// 	unsigned b = (me % num01);
	// 	MPI_Recv(public_key, public_key_len, MPI_CHAR, b, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	// }

	// printf("me1 == %d\n", me);
	// MPI_Barrier(MPI_COMM_WORLD);

	// MPI_Bcast(public_key, public_key_len, MPI_CHAR, 0, local_comm02);
	// printf("**************03 me = %d\n", me);
	for (size_t i = 0; i < 8; i++)
		hss_zeroize(&ctx[i], sizeof ctx[i]);

	// end = MPI_Wtime();
	// if(me == 0) printf("2 %lf\n", 1000 * (end - start));

	return true;
}

bool lm_ots_generate_public_key(
	param_set_t lm_ots_type,
	const unsigned char *I, /* Public key identifier */
	merkle_index_t q,       /* Diversification string, 4 bytes value */
	struct seed_derive *seed,
	unsigned char *public_key, size_t public_key_len)
{
	/* Look up the parameter set */
	unsigned h, n, w, p, ls;

	if (!lm_ots_look_up_parameter_set(lm_ots_type, &h, &n, &w, &p, &ls))
		return false;

	/* Start the hash that computes the final value */
	union hash_context public_ctx;

	hss_init_hash_context(h, &public_ctx);
	{
		unsigned char prehash_prefix[PBLC_PREFIX_LEN];
		memcpy(prehash_prefix + PBLC_I, I, I_LEN);
		put_bigendian(prehash_prefix + PBLC_Q, q, 4);
		SET_D(prehash_prefix + PBLC_D, D_PBLC);
		hss_update_hash_context(h, &public_ctx, prehash_prefix,
					PBLC_PREFIX_LEN);
	}

	int i, j;

	hss_seed_derive_set_j(seed, 0);
	int num = p / 8;
	int rem = p % 8;
	unsigned char buffer[(num + 1) * 8][ITER_MAX_LEN];
	union hash_context ctx[8];

	for (size_t i = 0; i < p; i++) {
		memcpy(buffer[i] + ITER_I, I, I_LEN);
		put_bigendian(buffer[i] + ITER_Q, q, 4);
		hss_seed_derive(buffer[i] + ITER_PREV, seed, i < p - 1);
		put_bigendian(buffer[i] + ITER_K, i, 2);
	}

	for (size_t i = 0; i < num; i++) {
		int offset = i * 8;
		for (j = 0; j < (1 << w) - 1; j++) {
			buffer[offset][ITER_J] = j;
			buffer[offset + 1][ITER_J] = j;
			buffer[offset + 2][ITER_J] = j;
			buffer[offset + 3][ITER_J] = j;
			buffer[offset + 4][ITER_J] = j;
			buffer[offset + 5][ITER_J] = j;
			buffer[offset + 6][ITER_J] = j;
			buffer[offset + 7][ITER_J] = j;
			hss_hash_ctx_eight(buffer[offset] + ITER_PREV, buffer[offset + 1] + ITER_PREV, buffer[offset + 2] + ITER_PREV, buffer[offset + 3] + ITER_PREV, buffer[offset + 4] + ITER_PREV, buffer[offset + 5] + ITER_PREV, buffer[offset + 6] + ITER_PREV, buffer[offset + 7] + ITER_PREV, h, buffer[offset], buffer[offset + 1], buffer[offset + 2], buffer[offset + 3], buffer[offset + 4], buffer[offset + 5], buffer[offset + 6], buffer[offset + 7], ITER_LEN(n));
		}
	}
	int m = p - rem;

	for (j = 0; j < (1 << w) - 1; j++) {
		for (size_t m = p - rem; m < p; m++)
			buffer[m][ITER_J] = j;
		m = p - rem;
		hss_hash_ctx_eight(buffer[m] + ITER_PREV, buffer[m + 1] + ITER_PREV, buffer[m + 2] + ITER_PREV, buffer[m + 3] + ITER_PREV, buffer[m + 4] + ITER_PREV, buffer[m + 5] + ITER_PREV, buffer[m + 6] + ITER_PREV, buffer[m + 7] + ITER_PREV, h, buffer[m], buffer[m + 1], buffer[m + 2], buffer[m + 3], buffer[m + 4], buffer[m + 5], buffer[m + 6], buffer[m + 7], ITER_LEN(n));
	}

	for (size_t i = 0; i < p; i++)
		hss_update_hash_context(h, &public_ctx, buffer[i] + ITER_PREV, n);
	/* And the result of the running hash is the public key */
	hss_finalize_hash_context(h, &public_ctx, public_key);

	for (size_t i = 0; i < 8; i++)
		hss_zeroize(&ctx[i], sizeof ctx[i]);
	return true;
}

#endif

#elif (defined(MPI_OTS))  //MPI--version
bool mpi_lm_ots_generate_public_key(
	param_set_t lm_ots_type,
	const unsigned char *I, /* Public key identifier */
	merkle_index_t q,       /* Diversification string, 4 bytes value */
	struct seed_derive *seed,
	unsigned char *public_key, size_t public_key_len)
{
	int me, proc_num;

	MPI_Comm_rank(MPI_COMM_WORLD, &me);             //进程号
	MPI_Comm_size(MPI_COMM_WORLD, &proc_num);       //进程数
	// if(me == 0) printf("1111111111\n");
	// double start, end;
	// if(me == 0) printf("22\n");
	// start = MPI_Wtime();
	unsigned int new_size;

	new_size = proc_num / global_level_nodes;

	// new_size = 2;
	MPI_Barrier(MPI_COMM_WORLD);

	// end = MPI_Wtime();
	// if(me == 0) printf("0 %lf\n", 1000 * (end - start));
	// start = MPI_Wtime();

	/* Look up the parameter set */
	unsigned h, n, w, p, ls;

	if (!lm_ots_look_up_parameter_set(lm_ots_type, &h, &n, &w, &p, &ls))
		return false;
	MPI_Barrier(MPI_COMM_WORLD);
	/* Start the hash that computes the final value */
	union hash_context public_ctx;

	hss_init_hash_context(h, &public_ctx);
	{
		unsigned char prehash_prefix[PBLC_PREFIX_LEN];
		memcpy(prehash_prefix + PBLC_I, I, I_LEN);
		put_bigendian(prehash_prefix + PBLC_Q, q, 4);
		SET_D(prehash_prefix + PBLC_D, D_PBLC);
		hss_update_hash_context(h, &public_ctx, prehash_prefix,
					PBLC_PREFIX_LEN);
	}
	int id;                                 //子通信域中新的进程号

#if (defined(Gather02)) || (defined(Igather02))
	MPI_Comm_rank(local_comm, &id);         //进程号
#else
	//假如没有划分通信子域
	id = me / global_level_nodes;
#endif
	// new_size是子通信域的总进程数
	int e = p / new_size;
	int f = p % new_size;
	unsigned int local = e + ((id < f) ? 1 : 0);
	unsigned int offset = id * e + ((id < f) ? id : f);

	hss_seed_derive_set_j(seed, 0);
	// int num = local / 8;
	// int rem = local % 8;
	// int m_num = 8;
	unsigned char buffer[p][ITER_MAX_LEN];
	union hash_context ctx;
	// unsigned char *bufs[m_num];
	// union hash_context ctx[8];
	// unsigned char empty[ITER_MAX_LEN];

	// for (size_t i = 0; i < 72; i++)
	// 	bufs[i] = &empty[0];

	for (size_t i = 0; i < p; i++) {
		memcpy(buffer[i] + ITER_I, I, I_LEN);
		put_bigendian(buffer[i] + ITER_Q, q, 4);
		hss_seed_derive(buffer[i] + ITER_PREV, seed, i < p - 1);
		put_bigendian(buffer[i] + ITER_K, i, 2);
	}
	for (size_t i = offset; i < offset + local; i++) {
		for (size_t j = 0; j < (1 << w) - 1; j++) {
			buffer[i][ITER_J] = j;

			hss_hash_ctx(buffer[i] + ITER_PREV, h, &ctx, buffer[i], ITER_LEN(n));
		}
	}
	//通信

	// end = MPI_Wtime();
	// if(me == 0) printf("1 %lf\n", 1000 * (end - start));
	// start = MPI_Wtime();
#if (defined(Gather02)) || (defined(Igather02))
	int r_counts[new_size];
	int r_dis[new_size];        //存起始位置

	for (int i = 0; i < f; i++) {
		r_counts[i] = (e + 1) * ITER_MAX_LEN;
		r_dis[i] = (i * e + i) * ITER_MAX_LEN;
	}
	for (int i = f; i < new_size; i++) {
		r_counts[i] = e * ITER_MAX_LEN;
		r_dis[i] = (i * e + f) * ITER_MAX_LEN;
	}
	MPI_Gatherv(buffer[offset], local * ITER_MAX_LEN, MPI_CHAR, buffer, r_counts, r_dis, MPI_CHAR, 0, local_comm);
#else /*思考：是否没有划通信子域，不可以用MPI_Gatherv？，选择使用send和recv*/
	// printf("************00me = %d\n", me);
	int a = me % global_level_nodes;
	if (id != 0)
		MPI_Send(buffer[offset], local * ITER_MAX_LEN, MPI_CHAR, a, 0, MPI_COMM_WORLD);
	if (id == 0) {
		if (f > 1) {
			for (size_t i = 1; i < f; i++) {
				unsigned int b = a + (i * global_level_nodes);
				MPI_Recv(buffer[i * e + i], (e + 1) * ITER_MAX_LEN, MPI_CHAR, b, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			}
		}
		for (size_t i = f; i < new_size; i++) {
			unsigned int b = a + (i * global_level_nodes);
			MPI_Recv(buffer[i * e + f], e * ITER_MAX_LEN, MPI_CHAR, b, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
	}

#endif

	for (size_t i = 0; i < p; i++)
		hss_update_hash_context(h, &public_ctx, buffer[i] + ITER_PREV, n);
	/* And the result of the running hash is the public key */
	hss_finalize_hash_context(h, &public_ctx, public_key);

	hss_zeroize(&ctx, sizeof ctx);
	// end = MPI_Wtime();
	// if(me == 0) printf("2 %lf\n", 1000 * (end - start));

	return true;
}
bool mpi_sign_lm_ots_generate_public_key(
	param_set_t lm_ots_type,
	const unsigned char *I, /* Public key identifier */
	merkle_index_t q,       /* Diversification string, 4 bytes value */
	struct seed_derive *seed,
	unsigned char *public_key, size_t public_key_len)
{
	int me, proc_num;

	MPI_Comm_rank(MPI_COMM_WORLD, &me);             //进程号
	MPI_Comm_size(MPI_COMM_WORLD, &proc_num);       //进程数
	// if(me == 0) printf("2222222222222\n");
	// if(me ==0) printf("11\n");
	// double start, end;

	// start = MPI_Wtime();

	// end = MPI_Wtime();
	// if(me == 0) printf("0 %lf\n", 1000 * (end - start));
	// start = MPI_Wtime();

	/* Look up the parameter set */
	unsigned h, n, w, p, ls;

	if (!lm_ots_look_up_parameter_set(lm_ots_type, &h, &n, &w, &p, &ls))
		return false;

	/* Start the hash that computes the final value */
	union hash_context public_ctx;

	hss_init_hash_context(h, &public_ctx);
	{
		unsigned char prehash_prefix[PBLC_PREFIX_LEN];
		memcpy(prehash_prefix + PBLC_I, I, I_LEN);
		put_bigendian(prehash_prefix + PBLC_Q, q, 4);
		SET_D(prehash_prefix + PBLC_D, D_PBLC);
		hss_update_hash_context(h, &public_ctx, prehash_prefix,
					PBLC_PREFIX_LEN);
	}
	int id, new_size;                  //子通信域中新的进程号

#ifdef Gather06

	MPI_Comm_rank(local_comm02, &id);               //进程号
	MPI_Comm_size(local_comm02, &new_size);         //进程数
#else
	id = me / num01;
	new_size = proc_num / num01;
#endif
	// new_size是子通信域的总进程数
	int e = p / new_size;
	int f = p % new_size;
	unsigned int local = e + ((id < f) ? 1 : 0);
	unsigned int offset = id * e + ((id < f) ? id : f);

	hss_seed_derive_set_j(seed, 0);
	// int num = local / 8;
	// int rem = local % 8;
	// int m_num = 8;
	unsigned char buffer[p][ITER_MAX_LEN];
	// unsigned char *bufs[m_num];
	union hash_context ctx;
	// union hash_context ctx[8];
	// unsigned char empty[ITER_MAX_LEN];

	// for (size_t i = 0; i < 72; i++)
	// 	bufs[i] = &empty[0];

	for (size_t i = 0; i < p; i++) {
		memcpy(buffer[i] + ITER_I, I, I_LEN);
		put_bigendian(buffer[i] + ITER_Q, q, 4);
		hss_seed_derive(buffer[i] + ITER_PREV, seed, i < p - 1);
		put_bigendian(buffer[i] + ITER_K, i, 2);
	}

	for (size_t i = offset; i < offset + local; i++) {
		for (size_t j = 0; j < (1 << w) - 1; j++) {
			buffer[i][ITER_J] = j;

			hss_hash_ctx(buffer[i] + ITER_PREV, h, &ctx, buffer[i], ITER_LEN(n));
		}
	}

#ifdef Gather06
	//通信
	int r_counts[new_size];
	int r_dis[new_size];        //存起始位置

	for (int i = 0; i < f; i++) {
		r_counts[i] = (e + 1) * ITER_MAX_LEN;
		r_dis[i] = (i * e + i) * ITER_MAX_LEN;
	}
	for (int i = f; i < new_size; i++) {
		r_counts[i] = e * ITER_MAX_LEN;
		r_dis[i] = (i * e + f) * ITER_MAX_LEN;
	}
	// end = MPI_Wtime();
	// if(me == 0) printf("1 %lf\n", 1000 * (end - start));
	// start = MPI_Wtime();

	MPI_Gatherv(buffer[offset], local * ITER_MAX_LEN, MPI_CHAR, buffer, r_counts, r_dis, MPI_CHAR, 0, local_comm02);
	if (id == 0) {
		for (size_t i = 0; i < p; i++)
			hss_update_hash_context(h, &public_ctx, buffer[i] + ITER_PREV, n);
		/* And the result of the running hash is the public key */
		hss_finalize_hash_context(h, &public_ctx, public_key);
	}
	MPI_Bcast(public_key, public_key_len, MPI_CHAR, 0, local_comm02);


#else

	int a = me % num01;
	if (id != 0)
		// printf("**************02 me = %d\n", me);
		MPI_Send(buffer[offset], local * ITER_MAX_LEN, MPI_CHAR, a, 0, MPI_COMM_WORLD);
	if (id == 0) {
		if (f > 1) {
			for (size_t i = 1; i < f; i++) {
				unsigned int b = a + (i * num01);
				MPI_Recv(buffer[i * e + i], (e + 1) * ITER_MAX_LEN, MPI_CHAR, b, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				// printf("**************03 me = %d\n", me);
			}
			for (size_t i = f; i < new_size; i++) {
				unsigned int b = a + (i * num01);
				MPI_Recv(buffer[i * e + f], e * ITER_MAX_LEN, MPI_CHAR, b, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				// printf("**************04 me = %d\n", me);
			}
		} else {
			for (size_t i = 1; i < new_size; i++) {
				unsigned int b = a + (i * num01);
				MPI_Recv(buffer[i * e + f], e * ITER_MAX_LEN, MPI_CHAR, b, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				// printf("**************04 me = %d\n", me);
			}
		}
	}
	if (id == 0) {
		for (size_t i = 0; i < p; i++)
			hss_update_hash_context(h, &public_ctx, buffer[i] + ITER_PREV, n);
		/* And the result of the running hash is the public key */
		hss_finalize_hash_context(h, &public_ctx, public_key);


		for (size_t i = 1; i < new_size; i++) {
			unsigned a = me + i * num01;
			MPI_Send(public_key, public_key_len, MPI_CHAR, a, 0, MPI_COMM_WORLD);
		}
	}
	if (id != 0) {
		unsigned b = (me % num01);
		MPI_Recv(public_key, public_key_len, MPI_CHAR, b, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}

#endif

	// MPI_Bcast(public_key, public_key_len, MPI_CHAR, 0, local_comm02);
	// printf("**************03 me = %d\n", me);

	hss_zeroize(&ctx, sizeof ctx);

	// end = MPI_Wtime();
	// if(me == 0) printf("2 %lf\n", 1000 * (end - start));

	return true;
}
bool lm_ots_generate_public_key(
	param_set_t lm_ots_type,
	const unsigned char *I, /* Public key identifier */
	merkle_index_t q,       /* Diversification string, 4 bytes value */
	struct seed_derive *seed,
	unsigned char *public_key, size_t public_key_len)
{
	//printf("111111111111111111111111111111111111\n");
	/* Look up the parameter set */
	unsigned h, n, w, p, ls;

	if (!lm_ots_look_up_parameter_set(lm_ots_type, &h, &n, &w, &p, &ls))
		return false;

	/* Start the hash that computes the final value */
	union hash_context public_ctx;

	hss_init_hash_context(h, &public_ctx);
	{
		unsigned char prehash_prefix[PBLC_PREFIX_LEN];
		memcpy(prehash_prefix + PBLC_I, I, I_LEN);
		put_bigendian(prehash_prefix + PBLC_Q, q, 4);
		SET_D(prehash_prefix + PBLC_D, D_PBLC);
		hss_update_hash_context(h, &public_ctx, prehash_prefix,
					PBLC_PREFIX_LEN);
	}

	/* Now generate the public key */
	/* This is where we spend the majority of the time during key gen and */
	/* signing operations; it would make sense to attempt to try to take */
	/* advantage of parallel (SIMD) hardware; even if we use it nowhere */
	/* else, we'd get a significant speed up */
	int i, j;

	unsigned char buf[ITER_MAX_LEN];

	memcpy(buf + ITER_I, I, I_LEN);
	put_bigendian(buf + ITER_Q, q, 4);
	union hash_context ctx;

	hss_seed_derive_set_j(seed, 0);

	for (i = 0; i < p; i++) {
		hss_seed_derive(buf + ITER_PREV, seed, i < p - 1);
		put_bigendian(buf + ITER_K, i, 2);
		/* We'll place j in the buffer below */
		for (j = 0; j < (1 << w) - 1; j++) {
			buf[ITER_J] = j;

			hss_hash_ctx(buf + ITER_PREV, h, &ctx, buf, ITER_LEN(n));
		}
		/* Include that in the hash */
		hss_update_hash_context(h, &public_ctx, buf + ITER_PREV, n);
	}

	/* And the result of the running hash is the public key */
	hss_finalize_hash_context(h, &public_ctx, public_key);

	hss_zeroize(&ctx, sizeof ctx);

	return true;
}
#else //original version
bool lm_ots_generate_public_key(
	param_set_t lm_ots_type,
	const unsigned char *I, /* Public key identifier */
	merkle_index_t q,       /* Diversification string, 4 bytes value */
	struct seed_derive *seed,
	unsigned char *public_key, size_t public_key_len)
{
	//printf("111111111111111111111111111111111111\n");
	/* Look up the parameter set */
	unsigned h, n, w, p, ls;

	if (!lm_ots_look_up_parameter_set(lm_ots_type, &h, &n, &w, &p, &ls))
		return false;

	/* Start the hash that computes the final value */
	union hash_context public_ctx;

	hss_init_hash_context(h, &public_ctx);
	{
		unsigned char prehash_prefix[PBLC_PREFIX_LEN];
		memcpy(prehash_prefix + PBLC_I, I, I_LEN);
		put_bigendian(prehash_prefix + PBLC_Q, q, 4);
		SET_D(prehash_prefix + PBLC_D, D_PBLC);
		hss_update_hash_context(h, &public_ctx, prehash_prefix,
					PBLC_PREFIX_LEN);
	}

	/* Now generate the public key */
	/* This is where we spend the majority of the time during key gen and */
	/* signing operations; it would make sense to attempt to try to take */
	/* advantage of parallel (SIMD) hardware; even if we use it nowhere */
	/* else, we'd get a significant speed up */
	int i, j;

	unsigned char buf[ITER_MAX_LEN];

	memcpy(buf + ITER_I, I, I_LEN);
	put_bigendian(buf + ITER_Q, q, 4);
	union hash_context ctx;

	hss_seed_derive_set_j(seed, 0);

	for (i = 0; i < p; i++) {
		hss_seed_derive(buf + ITER_PREV, seed, i < p - 1);
		put_bigendian(buf + ITER_K, i, 2);
		/* We'll place j in the buffer below */
		for (j = 0; j < (1 << w) - 1; j++) {
			buf[ITER_J] = j;

			hss_hash_ctx(buf + ITER_PREV, h, &ctx, buf, ITER_LEN(n));
		}
		/* Include that in the hash */
		hss_update_hash_context(h, &public_ctx, buf + ITER_PREV, n);
	}

	/* And the result of the running hash is the public key */
	hss_finalize_hash_context(h, &public_ctx, public_key);

	hss_zeroize(&ctx, sizeof ctx);

	return true;
}
#endif
/*
 * This generates the randomizer C.  We assume seed has been initialized to
 * the expected q value
 */
void lm_ots_generate_randomizer(unsigned char *c, unsigned n,
				struct seed_derive *seed)
{
	unsigned char randomizer[SEED_LEN];

	hss_seed_derive_set_j(seed, SEED_RANDOMIZER_INDEX);

	hss_seed_derive(randomizer, seed, false);

	memcpy(c, randomizer, n);
}

#if (defined(MPI_OTS))
double sign_result_OTS = 0;
double sign_count = 0;
bool lm_ots_generate_signature(
	param_set_t lm_ots_type,
	const unsigned char *I, /* Public key identifier */
	merkle_index_t q,       /* Diversification string, 4 bytes value */
	struct seed_derive *seed,
	const void *message, size_t message_len, bool prehashed,
	unsigned char *signature, size_t signature_len)
{
#if (defined(MPI_Time))
	double start, end;
	start = MPI_Wtime();
#endif
	int me, proc_num;

	MPI_Comm_rank(MPI_COMM_WORLD, &me);             //进程号
	MPI_Comm_size(MPI_COMM_WORLD, &proc_num);       //进程数

	/* Look up the parameter set */
	// struct timespec start, stop;
	//
	// //
	// clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start);
	unsigned h, n, w, p, ls;

	if (!lm_ots_look_up_parameter_set(lm_ots_type, &h, &n, &w, &p, &ls))
		return false;

	/* Check if we have enough room */
	if (signature_len < 4 + n + p * n) return false;

	/* Export the parameter set to the signature */
	put_bigendian(signature, lm_ots_type, 4);

	union hash_context ctx;

	/* Select the randomizer.  Note: we do this determanistically, because
	 * upper levels of the HSS tree sometimes sign the same message with the
	 * same index (between multiple reboots), hence we want to make sure that
	 * the randomizer for a particualr index is the same
	 * Also, if we're prehashed, we assume the caller has already selected it,
	 * and placed it into the siganture */

	if (!prehashed)
		lm_ots_generate_randomizer(signature + 4, n, seed);

	/* Compute the initial hash */
	unsigned char Q[MAX_HASH + 2];

	if (!prehashed) {
		hss_init_hash_context(h, &ctx);

		/* First, we hash the message prefix */
		/*将LMS密钥标识符I、LMS叶标识符q、值D\u MESG（0x8181）和随机化器C前置到消息*/
		unsigned char prefix[MESG_PREFIX_MAXLEN];
		memcpy(prefix + MESG_I, I, I_LEN);
		put_bigendian(prefix + MESG_Q, q, 4);
		SET_D(prefix + MESG_D, D_MESG);
		memcpy(prefix + MESG_C, signature + 4, n);
		hss_update_hash_context(h, &ctx, prefix, MESG_PREFIX_LEN(n));

		/* Then, the message */
		/*计算消息摘要*/
		hss_update_hash_context(h, &ctx, message, message_len);
		hss_finalize_hash_context(h, &ctx, Q);
	} else {
		memcpy(Q, message, n);
	}

	/* Append the checksum to the randomized hash */
	/*将哈希的校验和连接到哈希本身*/
	put_bigendian(&Q[n], lm_ots_compute_checksum(Q, n, w, ls), 2);

	hss_seed_derive_set_j(seed, 0);

	int steps[p];

	for (size_t j = 0; j < p; j++)
		steps[j] = lm_ots_coef(Q, j, w);

//方法二//

	int e = p / proc_num;
	int f = p % proc_num;
	unsigned int local = e + ((me < f) ? 1 : 0);
	unsigned int offset = me * e + ((me < f) ? me : f);

	// int m_num = 8;
	unsigned char buffer[p][ITER_MAX_LEN];

	// unsigned char *bufs[m_num];
	// unsigned char empty[ITER_MAX_LEN];

	// unsigned char tem_hash[n];
	//
	// hss_seed_derive(tem_hash, seed, 1);

	for (size_t i = 0; i < p; i++) {
		memcpy(buffer[i] + ITER_I, I, I_LEN);
		put_bigendian(buffer[i] + ITER_Q, q, 4);
		put_bigendian(buffer[i] + ITER_K, i, 2);
		hss_seed_derive(buffer[i] + ITER_PREV, seed, i < p - 1);
	}

	for (int i = offset; i < offset + local; i++) {
		unsigned a = lm_ots_coef(Q, i, w);
		unsigned j;
		for (j = 0; j < a; j++) {
			buffer[i][ITER_J] = j;
			hss_hash_ctx(buffer[i] + ITER_PREV, h, &ctx, buffer[i], ITER_LEN(n));
		}
	}
//通信
	int r_counts[proc_num];
	int r_dis[proc_num];//存起始位置

	for (int i = 0; i < f; i++) {
		r_counts[i] = (e + 1) * ITER_MAX_LEN;
		r_dis[i] = (i * e + i) * ITER_MAX_LEN;
	}
	for (int i = f; i < proc_num; i++) {
		r_counts[i] = e * ITER_MAX_LEN;
		r_dis[i] = (i * e + f) * ITER_MAX_LEN;
	}
	MPI_Gatherv(buffer[offset], local * ITER_MAX_LEN, MPI_CHAR, buffer, r_counts, r_dis, MPI_CHAR, 0, MPI_COMM_WORLD);
	// MPI_Allgatherv(buffer[offset], local * ITER_MAX_LEN, MPI_CHAR, buffer, r_counts, r_dis, MPI_CHAR, MPI_COMM_WORLD);

	for (int i = 0; i < p; i++)
		memcpy(&signature[4 + n + n * i], buffer[i] + ITER_PREV, n);

	hss_zeroize(&ctx, sizeof ctx);
	// clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &stop);
	// int me, proc_num;

	// MPI_Comm_rank(MPI_COMM_WORLD, &me);             //进程号
	// result += (stop.tv_sec - start.tv_sec) * 1e6 + (stop.tv_nsec - start.tv_nsec) / 1e3;
#if (defined(MPI_Time))
	end = MPI_Wtime();
	sign_count++;
	sign_result_OTS += (end - start);
	if (me == 0)
		printf("sign_count = %lf MPI sign_result_OTS = %.4lf msec %.4lf msec\n", sign_count, sign_result_OTS * 1e3, sign_result_OTS * 1e3 / sign_count);
#endif
	return true;
}

#elif (defined(three_parallel))
/*一次性签名这部分不区分256和512位宽*/
double sign_result_OTS = 0;
double sign_count = 0;
bool lm_ots_generate_signature(
	param_set_t lm_ots_type,
	const unsigned char *I, /* Public key identifier */
	merkle_index_t q,       /* Diversification string, 4 bytes value */
	struct seed_derive *seed,
	const void *message, size_t message_len, bool prehashed,
	unsigned char *signature, size_t signature_len)
{
#if (defined(MPI_Time))
	double start, end;
	start = MPI_Wtime();
#endif
	int me, proc_num;

	MPI_Comm_rank(MPI_COMM_WORLD, &me);             //进程号
	MPI_Comm_size(MPI_COMM_WORLD, &proc_num);       //进程数

	/* Look up the parameter set */
	// struct timespec start, stop;


	// clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start);
	unsigned h, n, w, p, ls;

	if (!lm_ots_look_up_parameter_set(lm_ots_type, &h, &n, &w, &p, &ls))
		return false;

	/* Check if we have enough room */
	if (signature_len < 4 + n + p * n) return false;

	/* Export the parameter set to the signature */
	put_bigendian(signature, lm_ots_type, 4);

	union hash_context ctx;

	/* Select the randomizer.  Note: we do this determanistically, because
	 * upper levels of the HSS tree sometimes sign the same message with the
	 * same index (between multiple reboots), hence we want to make sure that
	 * the randomizer for a particualr index is the same
	 * Also, if we're prehashed, we assume the caller has already selected it,
	 * and placed it into the siganture */

	if (!prehashed)
		lm_ots_generate_randomizer(signature + 4, n, seed);

	/* Compute the initial hash */
	unsigned char Q[MAX_HASH + 2];

	if (!prehashed) {
		hss_init_hash_context(h, &ctx);

		/* First, we hash the message prefix */
		/*将LMS密钥标识符I、LMS叶标识符q、值D\u MESG（0x8181）和随机化器C前置到消息*/
		unsigned char prefix[MESG_PREFIX_MAXLEN];
		memcpy(prefix + MESG_I, I, I_LEN);
		put_bigendian(prefix + MESG_Q, q, 4);
		SET_D(prefix + MESG_D, D_MESG);
		memcpy(prefix + MESG_C, signature + 4, n);
		hss_update_hash_context(h, &ctx, prefix, MESG_PREFIX_LEN(n));

		/* Then, the message */
		/*计算消息摘要*/
		hss_update_hash_context(h, &ctx, message, message_len);
		hss_finalize_hash_context(h, &ctx, Q);
	} else {
		memcpy(Q, message, n);
	}

	/* Append the checksum to the randomized hash */
	/*将哈希的校验和连接到哈希本身*/
	put_bigendian(&Q[n], lm_ots_compute_checksum(Q, n, w, ls), 2);

	hss_seed_derive_set_j(seed, 0);

	int steps[p];

	for (size_t j = 0; j < p; j++)
		steps[j] = lm_ots_coef(Q, j, w);

//方法二//

	int e = p / proc_num;
	int f = p % proc_num;
	unsigned int local = e + ((me < f) ? 1 : 0);
	unsigned int offset = me * e + ((me < f) ? me : f);

	int m_num = 8;
	unsigned char buffer[p][ITER_MAX_LEN];
	unsigned char *bufs[m_num];
	unsigned char empty[ITER_MAX_LEN];

	// unsigned char tem_hash[n];
	//
	// hss_seed_derive(tem_hash, seed, 1);

	for (size_t i = 0; i < p; i++) {
		memcpy(buffer[i] + ITER_I, I, I_LEN);
		put_bigendian(buffer[i] + ITER_Q, q, 4);
		put_bigendian(buffer[i] + ITER_K, i, 2);
		// hss_seed_derive(buffer[i] + ITER_PREV, seed, i < p - 1);
		// if (i < 2)
		// 	memcpy(buffer[i] + ITER_PREV, tem_hash, n);
		// else
		hss_seed_derive(buffer[i] + ITER_PREV, seed, i < p - 1);
	}

	int k[p];  //已经哈希次数

	memset(k, 0, p * sizeof(int));

	for (int i = 0; i < m_num; i++)
		bufs[i] = &empty[0];

	while (1) {
		int num = 0;            //0-7
		int target[m_num];      //下标
		for (size_t i = offset; i < offset + local; i++) {
			if (steps[i] > 0) {
				bufs[num] = buffer[i];
				target[num] = i;
				k[target[num]]++;
				num++;
			}
			if (num == m_num)
				break;
		}
		for (int i = 0; i < num; i++)
			bufs[i][ITER_J] = k[target[i]] - 1;

		hss_hash_ctx_eight(
			bufs[0] + ITER_PREV, bufs[1] + ITER_PREV, bufs[2] + ITER_PREV, bufs[3] + ITER_PREV,
			bufs[4] + ITER_PREV, bufs[5] + ITER_PREV, bufs[6] + ITER_PREV, bufs[7] + ITER_PREV,
			h,
			bufs[0], bufs[1], bufs[2], bufs[3], bufs[4], bufs[5], bufs[6], bufs[7], ITER_LEN(n));

		for (int i = 0; i < m_num; i++)
			bufs[i] = &empty[0];

		for (int i = 0; i < num; i++)
			steps[target[i]] -= 1;

		if (num == 0)
			break;
	}
#if (defined(Gather04))
	int r_counts[proc_num];
	int r_dis[proc_num];//存起始位置

	for (int i = 0; i < f; i++) {
		r_counts[i] = (e + 1) * ITER_MAX_LEN;
		r_dis[i] = (i * e + i) * ITER_MAX_LEN;
	}
	for (int i = f; i < proc_num; i++) {
		r_counts[i] = e * ITER_MAX_LEN;
		r_dis[i] = (i * e + f) * ITER_MAX_LEN;
	}
#if (defined(Communication_time))
	MPI_Barrier(MPI_COMM_WORLD);
	double t0, t1, t2, t_sum;
	t0 = MPI_Wtime();
#endif
	MPI_Gatherv(buffer[offset], local * ITER_MAX_LEN, MPI_CHAR, buffer, r_counts, r_dis, MPI_CHAR, 0, MPI_COMM_WORLD);
	// MPI_Allgatherv(buffer[offset], local * ITER_MAX_LEN, MPI_CHAR, buffer, r_counts, r_dis, MPI_CHAR, MPI_COMM_WORLD);
#if (defined(Communication_time))
	t1 = MPI_Wtime();
	t2 = t1 - t0;
	// printf("t2 = %f\n", t2 * 1000);

	MPI_Reduce(&t2, &t_sum, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	sum_time4 += (t_sum * 1000);
	if (me == 0)
		printf("LM-OTS_signature Gather_time = %.6lf***************\n", sum_time4);
#endif
#elif (defined(Igather04))
	int r_counts[proc_num];
	int r_dis[proc_num];//存起始位置

	for (int i = 0; i < f; i++) {
		r_counts[i] = (e + 1) * ITER_MAX_LEN;
		r_dis[i] = (i * e + i) * ITER_MAX_LEN;
	}
	for (int i = f; i < proc_num; i++) {
		r_counts[i] = e * ITER_MAX_LEN;
		r_dis[i] = (i * e + f) * ITER_MAX_LEN;
	}
#if (defined(Communication_time))
	MPI_Barrier(MPI_COMM_WORLD);
	double t0, t1, t2, t_sum;
	t0 = MPI_Wtime();
#endif
	MPI_Request re;
	// MPI_Iallgatherv(buffer[offset], local * ITER_MAX_LEN, MPI_CHAR, buffer, r_counts, r_dis, MPI_CHAR, MPI_COMM_WORLD, &re)
	MPI_Igatherv(buffer[offset], local * ITER_MAX_LEN, MPI_CHAR, buffer, r_counts, r_dis, MPI_CHAR, 0, MPI_COMM_WORLD, &re);
	MPI_Wait(&re, MPI_STATUS_IGNORE);
#if (defined(Communication_time))
	t1 = MPI_Wtime();
	t2 = t1 - t0;
	// printf("t2 = %f\n", t2 * 1000);

	MPI_Reduce(&t2, &t_sum, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	sum_time4 += (t_sum * 1000);
	if (me == 0)
		printf("LM-OTS_signature Igather_time = %.6lf***************\n", sum_time4);
#endif
#elif (defined(Send_Recv04))
/*点对点通信*/
#if (defined(Communication_time))
	MPI_Barrier(MPI_COMM_WORLD);
	double t0, t1, t2, t_sum;
	t0 = MPI_Wtime();
#endif
	if (me != 0 && proc_num > 1)
		MPI_Send(buffer[offset], local * ITER_MAX_LEN, MPI_CHAR, 0, 0, MPI_COMM_WORLD);

	if (me == 0 && proc_num > 1) {
		if (f > 1)
			for (size_t i = 1; i < f; i++)
				MPI_Recv(buffer[i * e + i], (e + 1) * ITER_MAX_LEN, MPI_CHAR, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		for (size_t i = f; i < proc_num; i++)
			MPI_Recv(buffer[i * e + f], e * ITER_MAX_LEN, MPI_CHAR, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
#if (defined(Communication_time))
	t1 = MPI_Wtime();
	t2 = t1 - t0;
	MPI_Reduce(&t2, &t_sum, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	sum_time4 += (t_sum * 1000);
	if (me == 0)
		printf("LM-OTS_signature Send_Recv_time = %.6lf***************\n", sum_time4);
#endif
#elif (defined(Isend_Irecv04))
/*非阻塞点对点通信*/
	MPI_Request request1, request2[proc_num];       //, request3[proc_num], request4;
#if (defined(Communication_time))
	MPI_Barrier(MPI_COMM_WORLD);                    // 同步
	double t0, t1, t2, t_sum;
	t0 = MPI_Wtime();
#endif
	if (me != 0 && proc_num > 1)
		MPI_Isend(buffer[offset], local * ITER_MAX_LEN, MPI_CHAR, 0, 0, MPI_COMM_WORLD, &request1);
	if (me == 0 && proc_num > 1) {
		if (f > 1)
			for (size_t i = 1; i < f; i++)
				MPI_Irecv(buffer[i * e + i], (e + 1) * ITER_MAX_LEN, MPI_CHAR, i, 0, MPI_COMM_WORLD, &request2[i]);
		for (size_t i = f; i < proc_num; i++)
			MPI_Irecv(buffer[i * e + f], e * ITER_MAX_LEN, MPI_CHAR, i, 0, MPI_COMM_WORLD, &request2[i]);
		// for (size_t i = 0; i < proc_num; i++) {
		// 	MPI_Isend(buffer, p * ITER_MAX_LEN, MPI_CHAR, i, 0, MPI_COMM_WORLD, &request3[i]);
		// }
	}
	// if (me != 0){
	// 	MPI_Irecv(buffer, p * ITER_MAX_LEN, MPI_CHAR, 0, 0, MPI_COMM_WORLD, &request4);
	// }
	if (me != 0)
		MPI_Waitall(1, &request1, MPI_STATUS_IGNORE);
	if (me == 0) {
		if (f > 1)
			MPI_Waitall(f - 1, &request2[1], MPI_STATUS_IGNORE);
		MPI_Waitall(proc_num - f, &request2[f], MPI_STATUS_IGNORE);
		// MPI_Waitall(proc_num, &request3, MPI_STATUS_IGNORE);
	}
	// if (me != 0){
	// 	MPI_Waitall(1, &request4, MPI_STATUS_IGNORE);
	// }

#if (defined(Communication_time))
	t1 = MPI_Wtime();
	t2 = t1 - t0;
	// printf("t2 = %f\n", t2 * 1000);

	MPI_Reduce(&t2, &t_sum, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	sum_time4 += (t_sum * 1000);
	if (me == 0)
		printf("************LM-OTS_signature Isend_Irecv Communication_time = %.6lf***************\n", sum_time4);
#endif
#endif
	for (int i = 0; i < p; i++)
		memcpy(&signature[4 + n + n * i], buffer[i] + ITER_PREV, n);

	hss_zeroize(&ctx, sizeof ctx);
	// clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &stop);
	// result += (stop.tv_sec - start.tv_sec) * 1e6 + (stop.tv_nsec - start.tv_nsec) / 1e3;
#if (defined(MPI_Time))
	end = MPI_Wtime();
	sign_count++;
	sign_result_OTS += (end - start);
	if (me == 0)
		printf("sign_count = %lf sign_result_OTS = %.4lf msec %.4lf msec\n", sign_count, sign_result_OTS * 1e3, sign_result_OTS * 1e3 / sign_count);
#endif
	return true;
}

#elif (defined(vectorization))
double sign_result_OTS = 0;
double sign_count = 0;

#if (defined(vector_512))

bool lm_ots_generate_signature(
	param_set_t lm_ots_type,
	const unsigned char *I, /* Public key identifier */
	merkle_index_t q,       /* Diversification string, 4 bytes value */
	struct seed_derive *seed,
	const void *message, size_t message_len, bool prehashed,
	unsigned char *signature, size_t signature_len)
{
	/* Look up the parameter set */
#if (defined(MPI_Time))
	double start, end;
	start = MPI_Wtime();
#endif
	int me, proc_num;

	MPI_Comm_rank(MPI_COMM_WORLD, &me);             //进程号
	MPI_Comm_size(MPI_COMM_WORLD, &proc_num);       //进程数
	// struct timespec start, stop;
	//
	//
	// clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start);
	unsigned h, n, w, p, ls;

	if (!lm_ots_look_up_parameter_set(lm_ots_type, &h, &n, &w, &p, &ls))
		return false;

	/* Check if we have enough room */
	if (signature_len < 4 + n + p * n) return false;

	/* Export the parameter set to the signature */
	put_bigendian(signature, lm_ots_type, 4);

	union hash_context ctx;


	if (!prehashed)
		lm_ots_generate_randomizer(signature + 4, n, seed);

	/* Compute the initial hash */
	unsigned char Q[MAX_HASH + 2];

	if (!prehashed) {
		hss_init_hash_context(h, &ctx);

		/* First, we hash the message prefix */
		/*将LMS密钥标识符I、LMS叶标识符q、值D\u MESG（0x8181）和随机化器C前置到消息*/
		unsigned char prefix[MESG_PREFIX_MAXLEN];
		memcpy(prefix + MESG_I, I, I_LEN);
		put_bigendian(prefix + MESG_Q, q, 4);
		SET_D(prefix + MESG_D, D_MESG);
		memcpy(prefix + MESG_C, signature + 4, n);
		hss_update_hash_context(h, &ctx, prefix, MESG_PREFIX_LEN(n));

		/* Then, the message */
		/*计算消息摘要*/
		hss_update_hash_context(h, &ctx, message, message_len);
		hss_finalize_hash_context(h, &ctx, Q);
	} else {
		memcpy(Q, message, n);
	}

	/* Append the checksum to the randomized hash */
	/*将哈希的校验和连接到哈希本身*/
	put_bigendian(&Q[n], lm_ots_compute_checksum(Q, n, w, ls), 2);

	hss_seed_derive_set_j(seed, 0);

	int steps[p];

	for (size_t j = 0; j < p; j++)
		steps[j] = lm_ots_coef(Q, j, w);


	int m_num = 16;
	unsigned char buffer[p][ITER_MAX_LEN];
	unsigned char *bufs[m_num];
	unsigned char empty[ITER_MAX_LEN];

	for (size_t i = 0; i < p; i++) {
		memcpy(buffer[i] + ITER_I, I, I_LEN);
		put_bigendian(buffer[i] + ITER_Q, q, 4);
		put_bigendian(buffer[i] + ITER_K, i, 2);
		hss_seed_derive(buffer[i] + ITER_PREV, seed, i < p - 1);
	}

	int k[p];  //已经哈希次数

	memset(k, 0, p * sizeof(int));
	while (1) {
		int num = 0;            //0-7
		int target[m_num];      //下标
		for (size_t i = 0; i < p; i++) {
			if (steps[i] > 0) {
				bufs[num] = buffer[i];
				target[num] = i;
				k[target[num]]++;
				num++;
			}
			if (num == m_num)
				break;
		}
		for (int i = 0; i < num; i++)
			bufs[i][ITER_J] = k[target[i]] - 1;

		hss_hash_ctx_sixteen(
			bufs[0] + ITER_PREV, bufs[1] + ITER_PREV, bufs[2] + ITER_PREV, bufs[3] + ITER_PREV,
			bufs[4] + ITER_PREV, bufs[5] + ITER_PREV, bufs[6] + ITER_PREV, bufs[7] + ITER_PREV,
			bufs[8] + ITER_PREV, bufs[9] + ITER_PREV, bufs[10] + ITER_PREV, bufs[11] + ITER_PREV,
			bufs[12] + ITER_PREV, bufs[13] + ITER_PREV, bufs[14] + ITER_PREV, bufs[15] + ITER_PREV,
			h,
			bufs[0], bufs[1], bufs[2], bufs[3], bufs[4], bufs[5], bufs[6], bufs[7],
			bufs[8], bufs[9], bufs[10], bufs[11], bufs[12], bufs[13], bufs[14], bufs[15],ITER_LEN(n));

		for (int i = 0; i < m_num; i++)
			bufs[i] = &empty[0];

		for (int i = 0; i < num; i++)
			steps[target[i]] -= 1;

		if (num == 0)
			break;
	}
	for (int i = 0; i < p; i++)
		memcpy(&signature[4 + n + n * i], buffer[i] + ITER_PREV, n);

/********************方法二*/
	hss_zeroize(&ctx, sizeof ctx);
	// clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &stop);
	// result += (stop.tv_sec - start.tv_sec) * 1e6 + (stop.tv_nsec - start.tv_nsec) / 1e3;
	// printf("Vector generate the VOTSsignature %.4lf msec\n", result / 1e3);
#if (defined(MPI_Time))
	end = MPI_Wtime();
	sign_count++;
	sign_result_OTS += (end - start);
	if (me == 0)
		printf("vector_512 sign_count = %lf Vector sign_result_OTS = %.4lf msec %.4lf msec\n", sign_count, sign_result_OTS * 1e3, sign_result_OTS * 1e3 / sign_count);
#endif
	// printf("^^^^^^^^\n");
	return true;
}


#else
bool lm_ots_generate_signature(
	param_set_t lm_ots_type,
	const unsigned char *I, /* Public key identifier */
	merkle_index_t q,       /* Diversification string, 4 bytes value */
	struct seed_derive *seed,
	const void *message, size_t message_len, bool prehashed,
	unsigned char *signature, size_t signature_len)
{
	/* Look up the parameter set */
#if (defined(MPI_Time))
	double start, end;
	start = MPI_Wtime();
#endif
	int me, proc_num;

	MPI_Comm_rank(MPI_COMM_WORLD, &me);             //进程号
	MPI_Comm_size(MPI_COMM_WORLD, &proc_num);       //进程数
	// struct timespec start, stop;
	//
	//
	// clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start);
	unsigned h, n, w, p, ls;

	if (!lm_ots_look_up_parameter_set(lm_ots_type, &h, &n, &w, &p, &ls))
		return false;

	/* Check if we have enough room */
	if (signature_len < 4 + n + p * n) return false;

	/* Export the parameter set to the signature */
	put_bigendian(signature, lm_ots_type, 4);

	union hash_context ctx;


	if (!prehashed)
		lm_ots_generate_randomizer(signature + 4, n, seed);

	/* Compute the initial hash */
	unsigned char Q[MAX_HASH + 2];

	if (!prehashed) {
		hss_init_hash_context(h, &ctx);

		/* First, we hash the message prefix */
		/*将LMS密钥标识符I、LMS叶标识符q、值D\u MESG（0x8181）和随机化器C前置到消息*/
		unsigned char prefix[MESG_PREFIX_MAXLEN];
		memcpy(prefix + MESG_I, I, I_LEN);
		put_bigendian(prefix + MESG_Q, q, 4);
		SET_D(prefix + MESG_D, D_MESG);
		memcpy(prefix + MESG_C, signature + 4, n);
		hss_update_hash_context(h, &ctx, prefix, MESG_PREFIX_LEN(n));

		/* Then, the message */
		/*计算消息摘要*/
		hss_update_hash_context(h, &ctx, message, message_len);
		hss_finalize_hash_context(h, &ctx, Q);
	} else {
		memcpy(Q, message, n);
	}

	/* Append the checksum to the randomized hash */
	/*将哈希的校验和连接到哈希本身*/
	put_bigendian(&Q[n], lm_ots_compute_checksum(Q, n, w, ls), 2);

	hss_seed_derive_set_j(seed, 0);

	int steps[p];

	for (size_t j = 0; j < p; j++)
		steps[j] = lm_ots_coef(Q, j, w);


//方法二//
	int m_num = 8;
	unsigned char buffer[p][ITER_MAX_LEN];
	unsigned char *bufs[m_num];
	unsigned char empty[ITER_MAX_LEN];

	for (size_t i = 0; i < p; i++) {
		memcpy(buffer[i] + ITER_I, I, I_LEN);
		put_bigendian(buffer[i] + ITER_Q, q, 4);
		put_bigendian(buffer[i] + ITER_K, i, 2);
		hss_seed_derive(buffer[i] + ITER_PREV, seed, i < p - 1);
	}

	int k[p];  //已经哈希次数

	memset(k, 0, p * sizeof(int));
	while (1) {
		int num = 0;            //0-7
		int target[m_num];      //下标
		for (size_t i = 0; i < p; i++) {
			if (steps[i] > 0) {
				bufs[num] = buffer[i];
				target[num] = i;
				k[target[num]]++;
				num++;
			}
			if (num == m_num)
				break;
		}
		for (int i = 0; i < num; i++)
			bufs[i][ITER_J] = k[target[i]] - 1;

		hss_hash_ctx_eight(
			bufs[0] + ITER_PREV, bufs[1] + ITER_PREV, bufs[2] + ITER_PREV, bufs[3] + ITER_PREV,
			bufs[4] + ITER_PREV, bufs[5] + ITER_PREV, bufs[6] + ITER_PREV, bufs[7] + ITER_PREV,
			h,
			bufs[0], bufs[1], bufs[2], bufs[3], bufs[4], bufs[5], bufs[6], bufs[7], ITER_LEN(n));

		for (int i = 0; i < m_num; i++)
			bufs[i] = &empty[0];

		for (int i = 0; i < num; i++)
			steps[target[i]] -= 1;

		if (num == 0)
			break;
	}
	for (int i = 0; i < p; i++)
		memcpy(&signature[4 + n + n * i], buffer[i] + ITER_PREV, n);

/********************方法二*/
	hss_zeroize(&ctx, sizeof ctx);
	// clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &stop);
	// result += (stop.tv_sec - start.tv_sec) * 1e6 + (stop.tv_nsec - start.tv_nsec) / 1e3;
	// printf("Vector generate the VOTSsignature %.4lf msec\n", result / 1e3);
#if (defined(MPI_Time))
	end = MPI_Wtime();
	sign_count++;
	sign_result_OTS += (end - start);
	if (me == 0)
		printf("Vector_256 sign_count = %lf Vector sign_result_OTS = %.4lf msec %.4lf msec\n", sign_count, sign_result_OTS * 1e3, sign_result_OTS * 1e3 / sign_count);
#endif
	// printf("^^^^^^^^\n");
	return true;
}
#endif


#else
double result = 0;
bool lm_ots_generate_signature(
	param_set_t lm_ots_type,
	const unsigned char *I, /* Public key identifier */
	merkle_index_t q,       /* Diversification string, 4 bytes value */
	struct seed_derive *seed,
	const void *message, size_t message_len, bool prehashed,
	unsigned char *signature, size_t signature_len)
{
	/* Look up the parameter set */
	struct timespec start, stop;


	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start);

	unsigned h, n, w, p, ls;

	if (!lm_ots_look_up_parameter_set(lm_ots_type, &h, &n, &w, &p, &ls))
		return false;

	/* Check if we have enough room */
	if (signature_len < 4 + n + p * n) return false;

	/* Export the parameter set to the signature */
	put_bigendian(signature, lm_ots_type, 4);

	union hash_context ctx;

	/* Select the randomizer.  Note: we do this determanistically, because
	 * upper levels of the HSS tree sometimes sign the same message with the
	 * same index (between multiple reboots), hence we want to make sure that
	 * the randomizer for a particualr index is the same
	 * Also, if we're prehashed, we assume the caller has already selected it,
	 * and placed it into the siganture */

	if (!prehashed)
		lm_ots_generate_randomizer(signature + 4, n, seed);

	/* Compute the initial hash */
	unsigned char Q[MAX_HASH + 2];

	if (!prehashed) {
		hss_init_hash_context(h, &ctx);

		/* First, we hash the message prefix */
		/*将LMS密钥标识符I、LMS叶标识符q、值D\u MESG（0x8181）和随机化器C前置到消息*/
		unsigned char prefix[MESG_PREFIX_MAXLEN];
		memcpy(prefix + MESG_I, I, I_LEN);
		put_bigendian(prefix + MESG_Q, q, 4);
		SET_D(prefix + MESG_D, D_MESG);
		memcpy(prefix + MESG_C, signature + 4, n);
		hss_update_hash_context(h, &ctx, prefix, MESG_PREFIX_LEN(n));

		/* Then, the message */
		/*计算消息摘要*/
		hss_update_hash_context(h, &ctx, message, message_len);
		hss_finalize_hash_context(h, &ctx, Q);
	} else {
		memcpy(Q, message, n);
	}

	/* Append the checksum to the randomized hash */
	/*将哈希的校验和连接到哈希本身*/
	put_bigendian(&Q[n], lm_ots_compute_checksum(Q, n, w, ls), 2);

	unsigned char tmp[ITER_MAX_LEN];

	/* Preset the parts of tmp that don't change */
	/*私钥tmp，将标识符I和q、密钥索引和散列计数器加入tmp，再进行散列a次得到OTS签名*/
	memcpy(tmp + ITER_I, I, I_LEN);
	put_bigendian(tmp + ITER_Q, q, 4);

	hss_seed_derive_set_j(seed, 0);

	for (int i = 0; i < p; i++) {
		put_bigendian(tmp + ITER_K, i, 2);
		hss_seed_derive(tmp + ITER_PREV, seed, i < p - 1);
		unsigned a = lm_ots_coef(Q, i, w);
		unsigned j;
		for (j = 0; j < a; j++) {
			tmp[ITER_J] = j;
			hss_hash_ctx(tmp + ITER_PREV, h, &ctx, tmp, ITER_LEN(n));
		}
		memcpy(&signature[4 + n + n * i], tmp + ITER_PREV, n);
	}

	hss_zeroize(&ctx, sizeof ctx);
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &stop);
	result += (stop.tv_sec - start.tv_sec) * 1e6 + (stop.tv_nsec - start.tv_nsec) / 1e3;
	// printf("original generate the OTSsignature %.4lf msec\n", result / 1e3);
	return true;
}
#endif
