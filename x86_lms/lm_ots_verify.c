/*
 * This is the code that implements the one-time-signature part of the LMS hash
 * based signatures
 */
#include <string.h>
#include "lm_ots_verify.h"
#include "lm_ots_common.h"
#include "hash.h"
#include "endian.h"
#include "common_defs.h"
#include "mpi.h"

/*
 * This validate a OTS signature for a message.  It doesn't actually use the
 * public key explicitly; instead, it just produces the root key, based on the
 * message; the caller is assumed to compare it to the expected value
 * Parameters:
 * - computed_public_key - where to place the reconstructed root.  It is
 *      assumed that the caller has allocated enough space
 * - I: the nonce value ("I") to use
 * - q: diversification string
 * - message - the message to verify
 * - message_len - the length of the message
 * - message_prehashed - true if the message has already undergone the initial
 *              (D_MESG) hash
 * - signature - the signature
 * - signature_len - the length of the signature
 * - parameter_set - what we expect the parameter set to be
 *
 * This returns true on successfully recomputing a root value; whether it is
 * the right one is something the caller would need to verify
 */
double sum_time5 = 0;

#if (defined(MPI_OTS))
bool lm_ots_validate_signature_compute(
	unsigned char *computed_public_key,
	const unsigned char *I, merkle_index_t q,
	const void *message, size_t message_len, bool message_prehashed,
	const unsigned char *signature, size_t signature_len,
	param_set_t expected_parameter_set)
{
// printf("1111111111\n");
	int me, proc_num;

	MPI_Comm_rank(MPI_COMM_WORLD, &me);             //进程号
	MPI_Comm_size(MPI_COMM_WORLD, &proc_num);       //进程数

	if (signature_len < 4) return false;            /* Ha, ha, very funny... */

	/* We don't trust the parameter set that's in the signature; verify it */
	param_set_t parameter_set = get_bigendian(signature, 4);

	if (parameter_set != expected_parameter_set)
		return false;

	unsigned h, n, w, p, ls;

	if (!lm_ots_look_up_parameter_set(parameter_set, &h, &n, &w, &p, &ls))
		return false;

	if (signature_len != 4 + n * (p + 1)) return false;

	const unsigned char *C = signature + 4;
	const unsigned char *y = C + n;

	unsigned char Q[MAX_HASH + 2];

	if (message_prehashed) {
		memcpy(Q, message, n);
	} else {
		union hash_context ctx;
		/* Compute the initial hash */
		hss_init_hash_context(h, &ctx);
		/* Hash the message prefix */
		{
			unsigned char prefix[MESG_PREFIX_MAXLEN];
			memcpy(prefix + MESG_I, I, I_LEN);
			put_bigendian(prefix + MESG_Q, q, 4);
			SET_D(prefix + MESG_D, D_MESG);
			memcpy(prefix + MESG_C, C, n);
			hss_update_hash_context(h, &ctx, prefix, MESG_PREFIX_LEN(n));
		}
		/* Then, the message */
		hss_update_hash_context(h, &ctx, message, message_len);

		hss_finalize_hash_context(h, &ctx, Q);
	}

	/* Append the checksum to the randomized hash */
	put_bigendian(&Q[n], lm_ots_compute_checksum(Q, n, w, ls), 2);

	/* And, start building the parts for the final hash */
	union hash_context final_ctx;

	hss_init_hash_context(h, &final_ctx);
	{
		unsigned char prehash_prefix[PBLC_PREFIX_LEN];
		memcpy(prehash_prefix + PBLC_I, I, I_LEN);
		put_bigendian(prehash_prefix + PBLC_Q, q, 4);
		SET_D(prehash_prefix + PBLC_D, D_PBLC);
		hss_update_hash_context(h, &final_ctx, prehash_prefix,
					PBLC_PREFIX_LEN);
	}

	unsigned max_digit = (1 << w) - 1;

	int e = p / proc_num;
	int f = p % proc_num;
	unsigned int local = e + ((me < f) ? 1 : 0);
	unsigned int offset = me * e + ((me < f) ? me : f);

	unsigned char buffer[p][ITER_MAX_LEN];


	for (size_t i = 0; i < p; i++) {
		memcpy(buffer[i] + ITER_I, I, I_LEN);
		put_bigendian(buffer[i] + ITER_Q, q, 4);
		put_bigendian(buffer[i] + ITER_K, i, 2);
		memcpy(buffer[i] + ITER_PREV, y + i * n, n);
	}
	for (size_t i = offset; i < offset + local; i++) {
		unsigned a = lm_ots_coef(Q, i, w);
		unsigned j;
		for (j = a; j < max_digit; j++) {
			union hash_context ctx;
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
	// MPI_Gatherv(buffer[offset], local * ITER_MAX_LEN, MPI_CHAR, buffer, r_counts, r_dis, MPI_CHAR, 0, MPI_COMM_WORLD);
	MPI_Allgatherv(buffer[offset], local * ITER_MAX_LEN, MPI_CHAR, buffer, r_counts, r_dis, MPI_CHAR, MPI_COMM_WORLD);


	for (int i = 0; i < p; i++)
		hss_update_hash_context(h, &final_ctx, buffer[i] + ITER_PREV, n);

	/*end向量化版本*/
	/* Ok, finalize the public key hash */
	hss_finalize_hash_context(h, &final_ctx, computed_public_key);
	// MPI_Bcast(computed_public_key, 32, MPI_CHAR, 0, MPI_COMM_WORLD);
	/*
	 * We succeeded in computing a root value; the caller will need to decide
	 * if the root we computed is actually the correct one
	 */
	return true;
}
#elif (defined(three_parallel))
#if (defined(vector_512))
bool lm_ots_validate_signature_compute(
	unsigned char *computed_public_key,
	const unsigned char *I, merkle_index_t q,
	const void *message, size_t message_len, bool message_prehashed,
	const unsigned char *signature, size_t signature_len,
	param_set_t expected_parameter_set)
{
	int me, proc_num;

	MPI_Comm_rank(MPI_COMM_WORLD, &me);             //进程号
	MPI_Comm_size(MPI_COMM_WORLD, &proc_num);       //进程数

	if (signature_len < 4) return false;            /* Ha, ha, very funny... */

	/* We don't trust the parameter set that's in the signature; verify it */
	param_set_t parameter_set = get_bigendian(signature, 4);

	if (parameter_set != expected_parameter_set)
		return false;

	unsigned h, n, w, p, ls;

	if (!lm_ots_look_up_parameter_set(parameter_set, &h, &n, &w, &p, &ls))
		return false;

	if (signature_len != 4 + n * (p + 1)) return false;

	const unsigned char *C = signature + 4;
	const unsigned char *y = C + n;

	unsigned char Q[MAX_HASH + 2];

	if (message_prehashed) {
		memcpy(Q, message, n);
	} else {
		union hash_context ctx;
		/* Compute the initial hash */
		hss_init_hash_context(h, &ctx);
		/* Hash the message prefix */
		{
			unsigned char prefix[MESG_PREFIX_MAXLEN];
			memcpy(prefix + MESG_I, I, I_LEN);
			put_bigendian(prefix + MESG_Q, q, 4);
			SET_D(prefix + MESG_D, D_MESG);
			memcpy(prefix + MESG_C, C, n);
			hss_update_hash_context(h, &ctx, prefix, MESG_PREFIX_LEN(n));
		}
		/* Then, the message */
		hss_update_hash_context(h, &ctx, message, message_len);

		hss_finalize_hash_context(h, &ctx, Q);
	}

	/* Append the checksum to the randomized hash */
	put_bigendian(&Q[n], lm_ots_compute_checksum(Q, n, w, ls), 2);

	/* And, start building the parts for the final hash */
	union hash_context final_ctx;

	hss_init_hash_context(h, &final_ctx);
	{
		unsigned char prehash_prefix[PBLC_PREFIX_LEN];
		memcpy(prehash_prefix + PBLC_I, I, I_LEN);
		put_bigendian(prehash_prefix + PBLC_Q, q, 4);
		SET_D(prehash_prefix + PBLC_D, D_PBLC);
		hss_update_hash_context(h, &final_ctx, prehash_prefix,
					PBLC_PREFIX_LEN);
	}

/*向量化版本*/
	int steps[p];
	int k[p];  //已经哈希次数
	unsigned max_digit = (1 << w) - 1;

	int e = p / proc_num;
	int f = p % proc_num;
	unsigned int local = e + ((me < f) ? 1 : 0);
	unsigned int offset = me * e + ((me < f) ? me : f);


	for (size_t j = 0; j < p; j++) {
		unsigned a = lm_ots_coef(Q, j, w);
		steps[j] = max_digit - a;
		k[j] = a;
	}

	int m_num = 16;
	unsigned char buffer[p][ITER_MAX_LEN];
	unsigned char *bufs[m_num];
	unsigned char empty[ITER_MAX_LEN];


	for (size_t i = 0; i < p; i++) {
		memcpy(buffer[i] + ITER_I, I, I_LEN);
		put_bigendian(buffer[i] + ITER_Q, q, 4);
		put_bigendian(buffer[i] + ITER_K, i, 2);
		memcpy(buffer[i] + ITER_PREV, y + i * n, n);
	}
	for (int i = 0; i < m_num; i++)
		bufs[i] = &empty[0];
	// if (me == 2) {
	// printf("%s\n", );
	while (1) {
		int num = 0;            //0-7
		int target[m_num];      //下标
		for (size_t i = offset; i < offset + local; i++) {
			// printf("steps[%d] = %d\n", i, steps[i]);
			if (steps[i] > 0) {
				bufs[num] = buffer[i];
				target[num] = i;
				k[target[num]]++;
				num++;
			}
			if (num == m_num)
				break;
		}
		// printf("num = %d\n",num );
		for (int i = 0; i < num; i++)
			bufs[i][ITER_J] = k[target[i]] - 1;

		hss_hash_ctx_sixteen(
			bufs[0] + ITER_PREV, bufs[1] + ITER_PREV, bufs[2] + ITER_PREV, bufs[3] + ITER_PREV,
			bufs[4] + ITER_PREV, bufs[5] + ITER_PREV, bufs[6] + ITER_PREV, bufs[7] + ITER_PREV,
			bufs[8] + ITER_PREV, bufs[9] + ITER_PREV, bufs[10] + ITER_PREV, bufs[11] + ITER_PREV,
			bufs[12] + ITER_PREV, bufs[13] + ITER_PREV, bufs[14] + ITER_PREV, bufs[15] + ITER_PREV,
			h,
			bufs[0], bufs[1], bufs[2], bufs[3], bufs[4], bufs[5], bufs[6], bufs[7],
			bufs[8], bufs[9], bufs[10], bufs[11], bufs[12], bufs[13], bufs[14], bufs[15],
			ITER_LEN(n));
		for (int i = 0; i < m_num; i++)
			bufs[i] = &empty[0];

		for (int i = 0; i < num; i++)
			steps[target[i]] -= 1;

		if (num == 0)
			break;
	}
	// }

	//
	// //通信
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
	// MPI_Gatherv(buffer[offset], local * ITER_MAX_LEN, MPI_CHAR, buffer, r_counts, r_dis, MPI_CHAR, 0, MPI_COMM_WORLD);
	MPI_Allgatherv(buffer[offset], local * ITER_MAX_LEN, MPI_CHAR, buffer, r_counts, r_dis, MPI_CHAR, MPI_COMM_WORLD);


	for (int i = 0; i < p; i++)
		hss_update_hash_context(h, &final_ctx, buffer[i] + ITER_PREV, n);

/*end向量化版本*/
	/* Ok, finalize the public key hash */
	hss_finalize_hash_context(h, &final_ctx, computed_public_key);
	// MPI_Bcast(computed_public_key, 32, MPI_CHAR, 0, MPI_COMM_WORLD);
	/*
	 * We succeeded in computing a root value; the caller will need to decide
	 * if the root we computed is actually the correct one
	 */
	return true;
}

#else
bool lm_ots_validate_signature_compute(
	unsigned char *computed_public_key,
	const unsigned char *I, merkle_index_t q,
	const void *message, size_t message_len, bool message_prehashed,
	const unsigned char *signature, size_t signature_len,
	param_set_t expected_parameter_set)
{
	int me, proc_num;

	MPI_Comm_rank(MPI_COMM_WORLD, &me);             //进程号
	MPI_Comm_size(MPI_COMM_WORLD, &proc_num);       //进程数

	if (signature_len < 4) return false;            /* Ha, ha, very funny... */

	/* We don't trust the parameter set that's in the signature; verify it */
	param_set_t parameter_set = get_bigendian(signature, 4);

	if (parameter_set != expected_parameter_set)
		return false;

	unsigned h, n, w, p, ls;

	if (!lm_ots_look_up_parameter_set(parameter_set, &h, &n, &w, &p, &ls))
		return false;

	if (signature_len != 4 + n * (p + 1)) return false;

	const unsigned char *C = signature + 4;
	const unsigned char *y = C + n;

	unsigned char Q[MAX_HASH + 2];

	if (message_prehashed) {
		memcpy(Q, message, n);
	} else {
		union hash_context ctx;
		/* Compute the initial hash */
		hss_init_hash_context(h, &ctx);
		/* Hash the message prefix */
		{
			unsigned char prefix[MESG_PREFIX_MAXLEN];
			memcpy(prefix + MESG_I, I, I_LEN);
			put_bigendian(prefix + MESG_Q, q, 4);
			SET_D(prefix + MESG_D, D_MESG);
			memcpy(prefix + MESG_C, C, n);
			hss_update_hash_context(h, &ctx, prefix, MESG_PREFIX_LEN(n));
		}
		/* Then, the message */
		hss_update_hash_context(h, &ctx, message, message_len);

		hss_finalize_hash_context(h, &ctx, Q);
	}

	/* Append the checksum to the randomized hash */
	put_bigendian(&Q[n], lm_ots_compute_checksum(Q, n, w, ls), 2);

	/* And, start building the parts for the final hash */
	union hash_context final_ctx;

	hss_init_hash_context(h, &final_ctx);
	{
		unsigned char prehash_prefix[PBLC_PREFIX_LEN];
		memcpy(prehash_prefix + PBLC_I, I, I_LEN);
		put_bigendian(prehash_prefix + PBLC_Q, q, 4);
		SET_D(prehash_prefix + PBLC_D, D_PBLC);
		hss_update_hash_context(h, &final_ctx, prehash_prefix,
					PBLC_PREFIX_LEN);
	}

/*向量化版本*/
	int steps[p];
	int k[p];  //已经哈希次数
	unsigned max_digit = (1 << w) - 1;

	int e = p / proc_num;
	int f = p % proc_num;
	unsigned int local = e + ((me < f) ? 1 : 0);
	unsigned int offset = me * e + ((me < f) ? me : f);


	for (size_t j = 0; j < p; j++) {
		unsigned a = lm_ots_coef(Q, j, w);
		steps[j] = max_digit - a;
		k[j] = a;
	}

	int m_num = 8;
	unsigned char buffer[p][ITER_MAX_LEN];
	unsigned char *bufs[m_num];
	unsigned char empty[ITER_MAX_LEN];


	for (size_t i = 0; i < p; i++) {
		memcpy(buffer[i] + ITER_I, I, I_LEN);
		put_bigendian(buffer[i] + ITER_Q, q, 4);
		put_bigendian(buffer[i] + ITER_K, i, 2);
		memcpy(buffer[i] + ITER_PREV, y + i * n, n);
	}
	for (int i = 0; i < m_num; i++)
		bufs[i] = &empty[0];
	// if (me == 2) {
	// printf("%s\n", );
	while (1) {
		int num = 0;            //0-7
		int target[m_num];      //下标
		for (size_t i = offset; i < offset + local; i++) {
			// printf("steps[%d] = %d\n", i, steps[i]);
			if (steps[i] > 0) {
				bufs[num] = buffer[i];
				target[num] = i;
				k[target[num]]++;
				num++;
			}
			if (num == m_num)
				break;
		}
		// printf("num = %d\n",num );
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
	// }

#if (defined(Gather05))
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
	// MPI_Gatherv(buffer[offset], local * ITER_MAX_LEN, MPI_CHAR, buffer, r_counts, r_dis, MPI_CHAR, 0, MPI_COMM_WORLD);
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
	sum_time5 += (t_sum * 1000);
	if (me == 0)
		printf("LM-OTS_verification Gather_time = %.6lf***************\n", sum_time5);
#endif
	for (int i = 0; i < p; i++)
		hss_update_hash_context(h, &final_ctx, buffer[i] + ITER_PREV, n);
	hss_finalize_hash_context(h, &final_ctx, computed_public_key);

#elif (defined(Igather05))
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
	// MPI_Gatherv(buffer[offset], local * ITER_MAX_LEN, MPI_CHAR, buffer, r_counts, r_dis, MPI_CHAR, 0, MPI_COMM_WORLD);
#if (defined(Communication_time))
	MPI_Barrier(MPI_COMM_WORLD);
	double t0, t1, t2, t_sum;
	t0 = MPI_Wtime();
#endif
	MPI_Request re;
	MPI_Igatherv(buffer[offset], local * ITER_MAX_LEN, MPI_CHAR, buffer, r_counts, r_dis, MPI_CHAR, 0, MPI_COMM_WORLD, &re);
	// MPI_Iallgatherv(buffer[offset], local * ITER_MAX_LEN, MPI_CHAR, buffer, r_counts, r_dis, MPI_CHAR, MPI_COMM_WORLD, &re);
	MPI_Wait(&re, MPI_STATUS_IGNORE);
#if (defined(Communication_time))
	t1 = MPI_Wtime();
	t2 = t1 - t0;
	// printf("t2 = %f\n", t2 * 1000);

	MPI_Reduce(&t2, &t_sum, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	sum_time5 += (t_sum * 1000);
	if (me == 0)
		printf("LM-OTS_verification Igather_time = %.6lf***************\n", sum_time5);
#endif
	for (int i = 0; i < p; i++)
		hss_update_hash_context(h, &final_ctx, buffer[i] + ITER_PREV, n);
	hss_finalize_hash_context(h, &final_ctx, computed_public_key);

#elif (defined(Send_Recv05))
#if (defined(Communication_time))
	MPI_Barrier(MPI_COMM_WORLD);
	double t0, t1, t2, t3, t4, t5, t_sum, t_sum1;
	t0 = MPI_Wtime();
#endif
	if (me != 0)
		MPI_Send(buffer[offset], local * ITER_MAX_LEN, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
	if (me == 0) {
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
	sum_time5 += (t_sum * 1000);
	if (me == 0)
		printf("LM-OTS_OTS_verification Send_Recv_time = %.6lf***************\n", sum_time5);
#endif
	if (me == 0) {
		for (int i = 0; i < p; i++)
			hss_update_hash_context(h, &final_ctx, buffer[i] + ITER_PREV, n);
		hss_finalize_hash_context(h, &final_ctx, computed_public_key);
	}

#elif (defined(Isend_Irecv05))
/*非阻塞点对点通信*/
	MPI_Request request1, request2[proc_num], request3[proc_num], request4;
#if (defined(Communication_time))
	MPI_Barrier(MPI_COMM_WORLD); // 同步
	double t0, t1, t2, t3, t4, t5, t_sum, t_sum1;
	t0 = MPI_Wtime();
#endif
	if (me != 0)
		MPI_Isend(buffer[offset], local * ITER_MAX_LEN, MPI_CHAR, 0, 0, MPI_COMM_WORLD, &request1);
	if (me == 0) {
		if (f > 1)
			for (size_t i = 1; i < f; i++)
				MPI_Irecv(buffer[i * e + i], (e + 1) * ITER_MAX_LEN, MPI_CHAR, i, 0, MPI_COMM_WORLD, &request2[i]);
		for (size_t i = f; i < proc_num; i++)
			MPI_Irecv(buffer[i * e + f], e * ITER_MAX_LEN, MPI_CHAR, i, 0, MPI_COMM_WORLD, &request2[i]);
	}
	if (me != 0)
		MPI_Waitall(1, &request1, MPI_STATUS_IGNORE);
	if (me == 0) {
		if (f > 1)
			MPI_Waitall(f - 1, &request2[1], MPI_STATUS_IGNORE);
		MPI_Waitall(proc_num - f, &request2[f], MPI_STATUS_IGNORE);
	}
#if (defined(Communication_time))
	t1 = MPI_Wtime();
	t2 = t1 - t0;
	// printf("t2 = %f\n", t2 * 1000);

	MPI_Reduce(&t2, &t_sum, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	sum_time5 += (t_sum * 1000);
	if (me == 0)
	    printf("************MAX LM-OTS_OTS_verification Isend_Irecv Communication_time = %.6lf***************\n", sum_time5);
#endif
	if (me == 0) {
		for (int i = 0; i < p; i++)
			hss_update_hash_context(h, &final_ctx, buffer[i] + ITER_PREV, n);
		hss_finalize_hash_context(h, &final_ctx, computed_public_key);
	}
#endif


/*end向量化版本*/
	return true;
}
#endif
#elif (defined(vectorization))
#if (defined(vector_512))
bool lm_ots_validate_signature_compute(
	unsigned char *computed_public_key,
	const unsigned char *I, merkle_index_t q,
	const void *message, size_t message_len, bool message_prehashed,
	const unsigned char *signature, size_t signature_len,
	param_set_t expected_parameter_set)
{

	if (signature_len < 4) return false; /* Ha, ha, very funny... */

	/* We don't trust the parameter set that's in the signature; verify it */
	param_set_t parameter_set = get_bigendian(signature, 4);

	if (parameter_set != expected_parameter_set)
		return false;

	unsigned h, n, w, p, ls;

	if (!lm_ots_look_up_parameter_set(parameter_set, &h, &n, &w, &p, &ls))
		return false;

	if (signature_len != 4 + n * (p + 1)) return false;

	const unsigned char *C = signature + 4;
	const unsigned char *y = C + n;

	unsigned char Q[MAX_HASH + 2];

	if (message_prehashed) {
		memcpy(Q, message, n);
	} else {
		union hash_context ctx;
		/* Compute the initial hash */
		hss_init_hash_context(h, &ctx);
		/* Hash the message prefix */
		{
			unsigned char prefix[MESG_PREFIX_MAXLEN];
			memcpy(prefix + MESG_I, I, I_LEN);
			put_bigendian(prefix + MESG_Q, q, 4);
			SET_D(prefix + MESG_D, D_MESG);
			memcpy(prefix + MESG_C, C, n);
			hss_update_hash_context(h, &ctx, prefix, MESG_PREFIX_LEN(n));
		}
		/* Then, the message */
		hss_update_hash_context(h, &ctx, message, message_len);

		hss_finalize_hash_context(h, &ctx, Q);
	}

	/* Append the checksum to the randomized hash */
	put_bigendian(&Q[n], lm_ots_compute_checksum(Q, n, w, ls), 2);

	/* And, start building the parts for the final hash */
	union hash_context final_ctx;

	hss_init_hash_context(h, &final_ctx);
	{
		unsigned char prehash_prefix[PBLC_PREFIX_LEN];
		memcpy(prehash_prefix + PBLC_I, I, I_LEN);
		put_bigendian(prehash_prefix + PBLC_Q, q, 4);
		SET_D(prehash_prefix + PBLC_D, D_PBLC);
		hss_update_hash_context(h, &final_ctx, prehash_prefix,
					PBLC_PREFIX_LEN);
	}

/*向量化版本*/
	int steps[p];
	int k[p]; //已经哈希次数
	unsigned max_digit = (1 << w) - 1;

	for (size_t j = 0; j < p; j++) {
		unsigned a = lm_ots_coef(Q, j, w);
		steps[j] = max_digit - a;
		k[j] = a;
	}

	int m_num = 16;
	unsigned char buffer[p][ITER_MAX_LEN];
	unsigned char *bufs[m_num];
	unsigned char empty[ITER_MAX_LEN];


	for (size_t i = 0; i < p; i++) {
		memcpy(buffer[i] + ITER_I, I, I_LEN);
		put_bigendian(buffer[i] + ITER_Q, q, 4);
		put_bigendian(buffer[i] + ITER_K, i, 2);
		memcpy(buffer[i] + ITER_PREV, y + i * n, n);
	}

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
		hss_update_hash_context(h, &final_ctx, buffer[i] + ITER_PREV, n);

/*end向量化版本*/
	/* Ok, finalize the public key hash */
	hss_finalize_hash_context(h, &final_ctx, computed_public_key);

	/*
	 * We succeeded in computing a root value; the caller will need to decide
	 * if the root we computed is actually the correct one
	 */
	return true;
}

#else
bool lm_ots_validate_signature_compute(
	unsigned char *computed_public_key,
	const unsigned char *I, merkle_index_t q,
	const void *message, size_t message_len, bool message_prehashed,
	const unsigned char *signature, size_t signature_len,
	param_set_t expected_parameter_set)
{

	if (signature_len < 4) return false; /* Ha, ha, very funny... */

	/* We don't trust the parameter set that's in the signature; verify it */
	param_set_t parameter_set = get_bigendian(signature, 4);

	if (parameter_set != expected_parameter_set)
		return false;

	unsigned h, n, w, p, ls;

	if (!lm_ots_look_up_parameter_set(parameter_set, &h, &n, &w, &p, &ls))
		return false;

	if (signature_len != 4 + n * (p + 1)) return false;

	const unsigned char *C = signature + 4;
	const unsigned char *y = C + n;

	unsigned char Q[MAX_HASH + 2];

	if (message_prehashed) {
		memcpy(Q, message, n);
	} else {
		union hash_context ctx;
		/* Compute the initial hash */
		hss_init_hash_context(h, &ctx);
		/* Hash the message prefix */
		{
			unsigned char prefix[MESG_PREFIX_MAXLEN];
			memcpy(prefix + MESG_I, I, I_LEN);
			put_bigendian(prefix + MESG_Q, q, 4);
			SET_D(prefix + MESG_D, D_MESG);
			memcpy(prefix + MESG_C, C, n);
			hss_update_hash_context(h, &ctx, prefix, MESG_PREFIX_LEN(n));
		}
		/* Then, the message */
		hss_update_hash_context(h, &ctx, message, message_len);

		hss_finalize_hash_context(h, &ctx, Q);
	}

	/* Append the checksum to the randomized hash */
	put_bigendian(&Q[n], lm_ots_compute_checksum(Q, n, w, ls), 2);

	/* And, start building the parts for the final hash */
	union hash_context final_ctx;

	hss_init_hash_context(h, &final_ctx);
	{
		unsigned char prehash_prefix[PBLC_PREFIX_LEN];
		memcpy(prehash_prefix + PBLC_I, I, I_LEN);
		put_bigendian(prehash_prefix + PBLC_Q, q, 4);
		SET_D(prehash_prefix + PBLC_D, D_PBLC);
		hss_update_hash_context(h, &final_ctx, prehash_prefix,
					PBLC_PREFIX_LEN);
	}

/*向量化版本*/
	int steps[p];
	int k[p]; //已经哈希次数
	unsigned max_digit = (1 << w) - 1;

	for (size_t j = 0; j < p; j++) {
		unsigned a = lm_ots_coef(Q, j, w);
		steps[j] = max_digit - a;
		k[j] = a;
	}

	int m_num = 8;
	unsigned char buffer[p][ITER_MAX_LEN];
	unsigned char *bufs[m_num];
	unsigned char empty[ITER_MAX_LEN];


	for (size_t i = 0; i < p; i++) {
		memcpy(buffer[i] + ITER_I, I, I_LEN);
		put_bigendian(buffer[i] + ITER_Q, q, 4);
		put_bigendian(buffer[i] + ITER_K, i, 2);
		memcpy(buffer[i] + ITER_PREV, y + i * n, n);
	}

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
		hss_update_hash_context(h, &final_ctx, buffer[i] + ITER_PREV, n);

/*end向量化版本*/
	/* Ok, finalize the public key hash */
	hss_finalize_hash_context(h, &final_ctx, computed_public_key);

	/*
	 * We succeeded in computing a root value; the caller will need to decide
	 * if the root we computed is actually the correct one
	 */
	return true;
}
#endif

#else
bool lm_ots_validate_signature_compute(
	unsigned char *computed_public_key,
	const unsigned char *I, merkle_index_t q,
	const void *message, size_t message_len, bool message_prehashed,
	const unsigned char *signature, size_t signature_len,
	param_set_t expected_parameter_set)
{
	if (signature_len < 4) return false; /* Ha, ha, very funny... */

	/* We don't trust the parameter set that's in the signature; verify it */
	param_set_t parameter_set = get_bigendian(signature, 4);

	if (parameter_set != expected_parameter_set)
		return false;

	unsigned h, n, w, p, ls;

	if (!lm_ots_look_up_parameter_set(parameter_set, &h, &n, &w, &p, &ls))
		return false;

	if (signature_len != 4 + n * (p + 1)) return false;

	const unsigned char *C = signature + 4;
	const unsigned char *y = C + n;

	unsigned char Q[MAX_HASH + 2];

	if (message_prehashed) {
		memcpy(Q, message, n);
	} else {
		union hash_context ctx;
		/* Compute the initial hash */
		hss_init_hash_context(h, &ctx);
		/* Hash the message prefix */
		{
			unsigned char prefix[MESG_PREFIX_MAXLEN];
			memcpy(prefix + MESG_I, I, I_LEN);
			put_bigendian(prefix + MESG_Q, q, 4);
			SET_D(prefix + MESG_D, D_MESG);
			memcpy(prefix + MESG_C, C, n);
			hss_update_hash_context(h, &ctx, prefix, MESG_PREFIX_LEN(n));
		}
		/* Then, the message */
		hss_update_hash_context(h, &ctx, message, message_len);

		hss_finalize_hash_context(h, &ctx, Q);
	}

	/* Append the checksum to the randomized hash */
	put_bigendian(&Q[n], lm_ots_compute_checksum(Q, n, w, ls), 2);

	/* And, start building the parts for the final hash */
	union hash_context final_ctx;

	hss_init_hash_context(h, &final_ctx);
	{
		unsigned char prehash_prefix[PBLC_PREFIX_LEN];
		memcpy(prehash_prefix + PBLC_I, I, I_LEN);
		put_bigendian(prehash_prefix + PBLC_Q, q, 4);
		SET_D(prehash_prefix + PBLC_D, D_PBLC);
		hss_update_hash_context(h, &final_ctx, prehash_prefix,
					PBLC_PREFIX_LEN);
	}


	int i;
	unsigned char tmp[ITER_MAX_LEN];

	/* Preset the parts of tmp that don't change */
	memcpy(tmp + ITER_I, I, I_LEN);
	put_bigendian(tmp + ITER_Q, q, 4);

	unsigned max_digit = (1 << w) - 1;

	for (i = 0; i < p; i++) {
		put_bigendian(tmp + ITER_K, i, 2);
		memcpy(tmp + ITER_PREV, y + i * n, n);
		unsigned a = lm_ots_coef(Q, i, w);
		unsigned j;
		for (j = a; j < max_digit; j++) {
			union hash_context ctx;
			tmp[ITER_J] = j;
			hss_hash_ctx(tmp + ITER_PREV, h, &ctx, tmp, ITER_LEN(n));
		}

		hss_update_hash_context(h, &final_ctx, tmp + ITER_PREV, n);
	}


	/* Ok, finalize the public key hash */
	hss_finalize_hash_context(h, &final_ctx, computed_public_key);

	/*
	 * We succeeded in computing a root value; the caller will need to decide
	 * if the root we computed is actually the correct one
	 */
	return true;
}
#endif
