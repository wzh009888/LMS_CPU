/*
 * This includes some computation methods that are shared between different
 * subsystems of the HSS signature package
 */
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "hss_internal.h"
#include "hss.h"
#include "hash.h"
#include "hss_thread.h"
#include "lm_ots_common.h"
#include "lm_ots.h"
#include "endian.h"
#include "hss_derive.h"
#include "mpi.h"
#include "time.h"


int num01;
int main_count = 0;
double sum_time3 = 0;
MPI_Comm local_comm02;
extern int global_level_nodes;
extern MPI_Comm local_comm001;
extern MPI_Comm local_comm04;
extern MPI_Comm local_comm08;
extern MPI_Comm local_comm16;
extern MPI_Comm local_comm64;


// MPI_Comm local_comm02;
/* Count the number of 1 bits at the end (lsbits) of the integer */
/* Do it in the obvious way; straightline code may be faster (no */
/* unpredictable jumps, which are costly), but that would be less scrutable */
/* (and this code is "fast enough") */
static int trailing_1_bits(merkle_index_t n)
{
	int i;

	for (i = 0; n& 1; n >>= 1, i++)
		;
	return i;
}

/*
 * Compute the value of an internal node within a Merkle tree
 */
double result_OTS = 0;
double sign_kg_OTS = 0;
double key_count = 0;
double sign_kg_count = 0;
double comm_count = 0;
#if (defined(three_parallel))
enum hss_error_code original_hss_compute_internal_node(unsigned char *		dest,
						       merkle_index_t		node_num,
						       const unsigned char *	seed,
						       param_set_t		lm_type,
						       param_set_t		lm_ots_type,
						       unsigned			h,
						       unsigned			leaf_level,
						       const unsigned char *	I)
{
	int me, proc_num;

	MPI_Comm_rank(MPI_COMM_WORLD, &me);             //进程号
	MPI_Comm_size(MPI_COMM_WORLD, &proc_num);       //进程数

	unsigned hash_size = hss_hash_length(h);

	/* We're store intermediate nodes here */
	unsigned char stack[MAX_HASH * MAX_MERKLE_HEIGHT];

	merkle_index_t tree_size = (merkle_index_t)1 << leaf_level;
	merkle_index_t r = node_num;
	unsigned levels_to_bottom = 0;

	if (r == 0) return hss_error_internal; /* So no to infinite loops */
	while (r < tree_size) {
		r <<= 1;
		levels_to_bottom++;
	}
	merkle_index_t q = r - tree_size;
	merkle_index_t i;
	unsigned ots_len = lm_ots_get_public_key_len(lm_ots_type);
	unsigned char pub_key[LEAF_MAX_LEN];

	memcpy(pub_key + LEAF_I, I, I_LEN);
	SET_D(pub_key + LEAF_D, D_LEAF);

	struct seed_derive derive;

	if (!hss_seed_derive_init(&derive, lm_type, lm_ots_type,
				  I, seed))
		return hss_error_bad_param_set;

	for (i = 0;; i++, r++, q++) {
		/* Generate the next OTS public key */
		hss_seed_derive_set_q(&derive, q);
#if (defined(MPI_Time)) || (defined(OTS_alltime))
		double start, end;
		start = MPI_Wtime();
#endif
		if (proc_num > global_level_nodes) {//进程数多于节点数
			// printf("0000000000000000\n" );
			if (!mpi_lm_ots_generate_public_key(lm_ots_type, I,
							    q, &derive, pub_key + LEAF_PK, ots_len))
				return hss_error_bad_param_set;         /* The only reason the above */
			/* could fail */
		} else {
			// printf("111111111111111111\n" );
			if (!lm_ots_generate_public_key(lm_ots_type, I,
							q, &derive, pub_key + LEAF_PK, ots_len))
				return hss_error_bad_param_set;         /* The only reason the above */
			/* could fail */
		}
#if (defined(MPI_Time)) || (defined(OTS_alltime))
		end = MPI_Wtime();
		key_count++;
		result_OTS += (end - start); //result单位是秒
#endif
		// if(me == 0) printf("1 %lf\n", 1000 * (end - start));
		// start = MPI_Wtime();

		/*
		 * For the subtree which this leaf node forms the final piece, put the
		 * destination to where we'll want it, either on the stack, or if this
		 * is the final piece, to where the caller specified
		 */
		unsigned char *current_buf;
		int stack_offset = trailing_1_bits(i);
		if (stack_offset == levels_to_bottom)
			current_buf = dest;
		else
			current_buf = &stack[stack_offset * hash_size];

		/* Hash it to form the leaf node */
		put_bigendian(pub_key + LEAF_R, r, 4);
		union hash_context ctx;
		hss_hash_ctx(current_buf, h, &ctx, pub_key, LEAF_LEN(hash_size));

		/* Work up the stack, combining right nodes with the left nodes */
		/* that we've already computed */

		unsigned sp;
		for (sp = 1; sp <= stack_offset; sp++) {
			hss_combine_internal_nodes(current_buf,
						   &stack[(sp - 1) * hash_size], current_buf,
						   h, I, hash_size,
						   r >> sp);
		}
		// end = MPI_Wtime();
		// if(me == 0) printf("2 %lf\n", 1000 * (end - start));

		/* We're not at a left branch, or at the target node */

		/* Because we've set current_buf to point to where we want to place */
		/* the result of this loop, we don't need to memcpy it */

		/* Check if this was the last leaf (and so we've just computed the */
		/* target node) */
		if (stack_offset == levels_to_bottom) {
			/* We're at the target node; the node we were asked to compute */
			/* We've already placed the value into dest, so we're all done */
			break;
		}
	}

	hss_seed_derive_done(&derive);

	return hss_error_none;
}

static enum hss_error_code hss_compute_internal_node(unsigned char *		dest,
						     merkle_index_t		node_num,
						     const unsigned char *	seed,
						     param_set_t		lm_type,
						     param_set_t		lm_ots_type,
						     unsigned			h,
						     unsigned			leaf_level,
						     const unsigned char *	I)
{
	int me, proc_num;

	MPI_Comm_rank(MPI_COMM_WORLD, &me);             //进程号
	MPI_Comm_size(MPI_COMM_WORLD, &proc_num);       //进程数

	main_count++;
	// struct timespec start, stop;

	// clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start);

	unsigned hash_size = hss_hash_length(h);

	/* We're store intermediate nodes here */
	unsigned char stack[MAX_HASH * MAX_MERKLE_HEIGHT];

	merkle_index_t tree_size = (merkle_index_t)1 << leaf_level;

	merkle_index_t r = node_num;
	unsigned levels_to_bottom = 0;
	//每个进程计算的树高和叶节点数
	unsigned new_levels_to_bottom = 0;
	int leaf_node = 0;

	if (r == 0) return hss_error_internal; /* So no to infinite loops */
	while (r < tree_size) {
		r <<= 1;
		levels_to_bottom++;
	}
	// if(me == 0) printf("me = %d r = %u levels_to_bottom = %u\n",me, r, levels_to_bottom );
	merkle_index_t q = r - tree_size;
	merkle_index_t i;
	unsigned ots_len = lm_ots_get_public_key_len(lm_ots_type);
	unsigned char pub_key[LEAF_MAX_LEN];

	memcpy(pub_key + LEAF_I, I, I_LEN);
	SET_D(pub_key + LEAF_D, D_LEAF);

	struct seed_derive derive;

	if (!hss_seed_derive_init(&derive, lm_type, lm_ots_type,
				  I, seed))
		return hss_error_bad_param_set;

	/*MPI*/

	unsigned char sum_current_buf[proc_num][32];

	num01 = pow(2, levels_to_bottom);
	int max_proc_num = num01;

	// if (me == 0) printf("num01 = %d\n", num01);
#ifdef Gather06
	if (proc_num > num01) {
		switch (num01) {
		case 1: MPI_Comm_dup(local_comm001, &local_comm02);
			break;
		case 4: MPI_Comm_dup(local_comm04, &local_comm02);
			break;
		case 8: MPI_Comm_dup(local_comm08, &local_comm02);
			break;
		case 16: MPI_Comm_dup(local_comm16, &local_comm02);
			break;
		case 64: MPI_Comm_dup(local_comm64, &local_comm02);
			break;
		default: break;
		}
	}
#endif
	int id;

	// if ((proc_num % 2) == 0){ //进程数为2的幂次情况
	if (proc_num > num01) {
		id = me % num01;
		leaf_node = 1;
		// printf("me = %d id = %d \n", me, id);
	} else {
		id = me;
		leaf_node = num01 / proc_num;
	}
	// }
	// else{	//进程数不是2的幂次情况
	// 	if (proc_num > num01) { //
	// 		unsigned int e =proc_num / num01; //1
	// 		unsigned int f = proc_num % num01; //192
	// 		unsigned int f_pro = e + 1;	//前f个节点每个节点e+1个进程算 ,剩下的num01 - f 个节点 每个节点e个进程算
	// 		unsigned int multipro = f * (e + 1);
	// 		if (me < multipro){ //前384个进程每两个一组算192个节点，有奇数问题吗？？？ 应该没有，只有一个进程算多个节点的情况有奇数节点数问题
	// 			id = me % f; //id值指示me进程算第id个叶节点
	// 		}
	// 		else{ //后256 - 192 = 64 个节点 每个节点分1个进程
	// 			id = (me - multipro) % (num01 - f) + f;
	// 		}
	// 		leaf_node = 1;
	// 		// printf("me = %d id = %d \n", me, id);
	// 	} else { //一个进程算多个节点
	// 		id = me; //id值指示me进程算第id个分支节点
	// 		// unsigned int min_2 = proc_num - proc_num % 2;
	// 		unsigned int e = (num01 / 2) / proc_num;
	// 		unsigned int f = (num01 / 2) % proc_num;
	// 		unsigned int local = e + ((id < f) ? 1 : 0);
	// 		leaf_node = local * 2;
	// 	}
	// }

	new_levels_to_bottom = log(leaf_node) / log(2);
	r = r + id * leaf_node;
	q = q + id * leaf_node;
	// printf("(((((((((((((())))))))))))))\n");
	for (i = 0;; i++, r++, q++) {
		/* Generate the next OTS public key */

		hss_seed_derive_set_q(&derive, q);
#if (defined(OTS_alltime))
		double start, end;
		start = MPI_Wtime();
#endif
		if (proc_num > num01) {      //进程数多于节点数
			if (!mpi_sign_lm_ots_generate_public_key(lm_ots_type, I,
								 q, &derive, pub_key + LEAF_PK, ots_len))
				return hss_error_bad_param_set;                 /* The only reason the above */
			/* could fail */
		} else {
			if (!lm_ots_generate_public_key(lm_ots_type, I,
							q, &derive, pub_key + LEAF_PK, ots_len))
				return hss_error_bad_param_set;                 /* The only reason the above */
			/* could fail */
		}
#if (defined(OTS_alltime))
		end = MPI_Wtime();
		sign_kg_count++;
		sign_kg_OTS += (end - start);
#endif
		unsigned char *current_buf;
		int stack_offset = trailing_1_bits(i);
		if (stack_offset == new_levels_to_bottom)
			current_buf = dest;
		else
			current_buf = &stack[stack_offset * hash_size];

		put_bigendian(pub_key + LEAF_R, r, 4);
		union hash_context ctx;
		hss_hash_ctx(current_buf, h, &ctx, pub_key, LEAF_LEN(hash_size));

		unsigned sp;
		for (sp = 1; sp <= stack_offset; sp++) {
			hss_combine_internal_nodes(current_buf,
						   &stack[(sp - 1) * hash_size], current_buf,
						   h, I, hash_size,
						   r >> sp);
		}

		if (stack_offset == new_levels_to_bottom)
			break;
	}

	// printf("!!!!!!!!!!me = %d\n", me);
	if (levels_to_bottom != 0) { //需通信后由0号进程接着构造子树
		// 划分通信域
		max_proc_num = max_proc_num > proc_num ? proc_num : max_proc_num;
#if (defined(Gather03)) || (defined(Igather03))
		MPI_Comm local_comm9;
		if (max_proc_num != proc_num) { //进程数多于num01,会有多个进程构建一个叶节点，故划分通信子域
			int local_rank[max_proc_num];
			MPI_Group World_Group, local_group;
			MPI_Comm_group(MPI_COMM_WORLD, &World_Group);

			for (int i = 0; i < max_proc_num; i++)
				local_rank[i] = i;

			MPI_Group_incl(World_Group, max_proc_num, local_rank, &local_group);
			MPI_Comm_create(MPI_COMM_WORLD, local_group, &local_comm9);
		}

#endif
#if (defined(Gather03))
#if (defined(Communication_time))
		MPI_Barrier(MPI_COMM_WORLD); // 同步
		double t0, t1, t2, t_sum;
		t0 = MPI_Wtime();
#endif
		if (max_proc_num != proc_num) { //进程数多于num01,会有多个进程构建一个叶节点，故划分通信子域
			if (me < max_proc_num)
				MPI_Gather(dest, 32, MPI_CHAR, sum_current_buf[0], 32, MPI_CHAR, 0, local_comm9);
		} else { //进程数少于或等于num01，不存在三级并行，直接在全局通信域中Gather即可
			// MPI_Request re;
			MPI_Gather(dest, 32, MPI_CHAR, sum_current_buf[0], 32, MPI_CHAR, 0, MPI_COMM_WORLD);
			// MPI_Wait(&re, MPI_STATUS_IGNORE);
		}
#if (defined(Communication_time))
		t1 = MPI_Wtime();
		t2 = t1 - t0;
		// printf("t2 = %f\n", t2 * 1000);

		MPI_Reduce(&t2, &t_sum, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
		if (me == 0) {
			sum_time3 += (t_sum * 1000);
			comm_count++;
			printf("************sum Gather Communication_time = %.6lf  %f ***************\n", sum_time3, comm_count);
		}
#endif
#elif (defined(Igather03))
#if (defined(Communication_time))
		MPI_Barrier(MPI_COMM_WORLD); // 同步
		double t0, t1, t2, t_sum;
		t0 = MPI_Wtime();
#endif
		if (max_proc_num != proc_num) { //进程数多于num01,会有多个进程构建一个叶节点，故划分通信子域
			if (me < max_proc_num) {
				MPI_Request re;
				MPI_Igather(dest, 32, MPI_CHAR, sum_current_buf[0], 32, MPI_CHAR, 0, local_comm9, &re);
			}
		} else { //进程数少于或等于num01，不存在三级并行，直接在全局通信域中Gather即可
			MPI_Request re;
			MPI_Igather(dest, 32, MPI_CHAR, sum_current_buf[0], 32, MPI_CHAR, 0, MPI_COMM_WORLD, &re);
			MPI_Wait(&re, MPI_STATUS_IGNORE);
		}
#if (defined(Communication_time))
		t1 = MPI_Wtime();
		t2 = t1 - t0;
		// printf("t2 = %f\n", t2 * 1000);

		MPI_Reduce(&t2, &t_sum, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
		if (me == 0) {
			sum_time3 += (t_sum * 1000);
			comm_count++;
			printf("************sum Igather Communication_time = %.6lf  %f ***************\n", sum_time3, comm_count);
		}
#endif
#elif (defined(Send_Recv03))
#if (defined(Communication_time))
		MPI_Barrier(MPI_COMM_WORLD); // 同步
		double t0, t1, t2, t_sum;
		t0 = MPI_Wtime();
#endif
		if (me < num01 && me != 0)
			MPI_Send(dest, 32, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
#if (defined(Communication_time))
		t1 = MPI_Wtime();
		t2 = t1 - t0;
#endif
#elif (defined(Isend_Irecv03))
/*非阻塞通信*/
#if (defined(Communication_time))
		MPI_Barrier(MPI_COMM_WORLD); // 同步
		double t0, t1, t2, t_sum;
		t0 = MPI_Wtime();
#endif
		MPI_Request request[max_proc_num - 1];
		if (me < num01 && me != 0)
			MPI_Send(dest, 32, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
		if (me == 0) {
			for (size_t i = 1; i < max_proc_num; i++)
				MPI_Irecv(sum_current_buf[i], 32, MPI_CHAR, i, 0, MPI_COMM_WORLD, &request[i - 1]);
			// MPI_Waitall(max_proc_num - 1, request, MPI_STATUS_IGNORE);
		}
#if (defined(Communication_time))
		t1 = MPI_Wtime();
		t2 = t1 - t0;
#endif
#endif
		// printf("&&&&&&&&&&me = %d\n", me);
		if (me == 0) {
			// main_count++;
			// MPI_Request request[max_proc_num - 1];
			memcpy(sum_current_buf[0], dest, 32);
			memcpy(dest, "0", 32);

			unsigned levels_to_bottom_0 = levels_to_bottom - new_levels_to_bottom;
			r = r >> new_levels_to_bottom;
			for (i = 0;; i++, r++) {
#if (defined(Send_Recv03))
#if (defined(Communication_time))
				t0 = MPI_Wtime();
#endif
				if (i != 0)
					MPI_Recv(sum_current_buf[i], 32, MPI_CHAR, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
#if (defined(Communication_time))
				t1 = MPI_Wtime();
				t2 = t1 - t0;
#endif
#elif (defined(Isend_Irecv03))
#if (defined(Communication_time))
				t0 = MPI_Wtime();
#endif
				if (i != 0)
					MPI_Wait(&request[i - 1], MPI_STATUS_IGNORE);
#if (defined(Communication_time))
				t1 = MPI_Wtime();
				t2 += (t1 - t0);
#endif
#endif
				/* Generate the next OTS public key */

				unsigned char *new_current_buf;
				int stack_offset = trailing_1_bits(i);
				if (stack_offset == levels_to_bottom_0)
					new_current_buf = dest;
				else
					new_current_buf = &stack[stack_offset * hash_size];

				memcpy(new_current_buf, sum_current_buf[i], 32);

				unsigned sp;
				for (sp = 1; sp <= stack_offset; sp++) {
					hss_combine_internal_nodes(new_current_buf,
								   &stack[(sp - 1) * hash_size], new_current_buf,
								   h, I, hash_size,
								   r >> sp);
				}

				if (stack_offset == levels_to_bottom_0)
					break;

				// printf("^^^^^^^^^^^^^^\n");
			}
		}
		// printf("8888888888888 me = %d\n", me);
#if (defined(Send_Recv03))
#if (defined(Communication_time))
		MPI_Reduce(&t2, &t_sum, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
		sum_time3 += (t_sum * 1000);
		if (me == 0) {
			comm_count++;
			printf("************MAX send_recv_time = %.6lf  %f***************\n", sum_time3, comm_count);
		}
#endif
#endif
#if (defined(Isend_Irecv03))
#if (defined(Communication_time))
		MPI_Reduce(&t2, &t_sum, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
		sum_time3 += (t_sum * 1000);
		if (me == 0) {
			comm_count++;
			printf("************MAX Isend_Irecv_time = %.6lf  %f***************\n", sum_time3, comm_count);
		}
#endif
#endif
		MPI_Bcast(dest, 32, MPI_CHAR, 0, MPI_COMM_WORLD);
		// printf("11111111111111111 me %d \n", me);
	}
	// MPI_Barrier(MPI_COMM_WORLD); // 同步

	hss_seed_derive_done(&derive);

	// clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &stop);
	// result += (stop.tv_sec - start.tv_sec) * 1e6 + (stop.tv_nsec - start.tv_nsec) / 1e3;
	// if (me == 0) printf("generate the nodee %.2lf msec\n", result / 1e3);
	// printf("+++++++++++++++++++me = %d \n", me);
	return hss_error_none;
}
// #endif
#elif (defined(parallel_tree))
enum hss_error_code original_hss_compute_internal_node(unsigned char *		dest,
						       merkle_index_t		node_num,
						       const unsigned char *	seed,
						       param_set_t		lm_type,
						       param_set_t		lm_ots_type,
						       unsigned			h,
						       unsigned			leaf_level,
						       const unsigned char *	I)
{
	unsigned hash_size = hss_hash_length(h);

	/* We're store intermediate nodes here */
	unsigned char stack[MAX_HASH * MAX_MERKLE_HEIGHT];

	merkle_index_t tree_size = (merkle_index_t)1 << leaf_level;
	merkle_index_t r = node_num;
	unsigned levels_to_bottom = 0;

	if (r == 0) return hss_error_internal; /* So no to infinite loops */
	while (r < tree_size) {
		r <<= 1;
		levels_to_bottom++;
	}
	merkle_index_t q = r - tree_size;

	merkle_index_t i;
	unsigned ots_len = lm_ots_get_public_key_len(lm_ots_type);
	unsigned char pub_key[LEAF_MAX_LEN];

	memcpy(pub_key + LEAF_I, I, I_LEN);
	SET_D(pub_key + LEAF_D, D_LEAF);

	struct seed_derive derive;

	if (!hss_seed_derive_init(&derive, lm_type, lm_ots_type,
				  I, seed))
		return hss_error_bad_param_set;

	for (i = 0;; i++, r++, q++) {
		/* Generate the next OTS public key */
		hss_seed_derive_set_q(&derive, q);
#if (defined(MPI_Time)) || (defined(OTS_alltime))
		double start, end;
		start = MPI_Wtime();
#endif
		if (!lm_ots_generate_public_key(lm_ots_type, I,
						q, &derive, pub_key + LEAF_PK, ots_len))
			return hss_error_bad_param_set;         /* The only reason the above */
		                                                /* could fail */

#if (defined(MPI_Time)) || (defined(OTS_alltime))
		end = MPI_Wtime();
		key_count++;
		result_OTS += (end - start);
#endif
		/*
		 * For the subtree which this leaf node forms the final piece, put the
		 * destination to where we'll want it, either on the stack, or if this
		 * is the final piece, to where the caller specified
		 */
		unsigned char *current_buf;
		int stack_offset = trailing_1_bits(i);
		if (stack_offset == levels_to_bottom)
			current_buf = dest;
		else
			current_buf = &stack[stack_offset * hash_size];

		/* Hash it to form the leaf node */
		put_bigendian(pub_key + LEAF_R, r, 4);
		union hash_context ctx;
		hss_hash_ctx(current_buf, h, &ctx, pub_key, LEAF_LEN(hash_size));

		/* Work up the stack, combining right nodes with the left nodes */
		/* that we've already computed */
		unsigned sp;
		for (sp = 1; sp <= stack_offset; sp++) {
			hss_combine_internal_nodes(current_buf,
						   &stack[(sp - 1) * hash_size], current_buf,
						   h, I, hash_size,
						   r >> sp);
		}

		/* We're not at a left branch, or at the target node */

		/* Because we've set current_buf to point to where we want to place */
		/* the result of this loop, we don't need to memcpy it */

		/* Check if this was the last leaf (and so we've just computed the */
		/* target node) */
		if (stack_offset == levels_to_bottom) {
			/* We're at the target node; the node we were asked to compute */
			/* We've already placed the value into dest, so we're all done */
			break;
		}
	}

	hss_seed_derive_done(&derive);

	return hss_error_none;
}

static enum hss_error_code hss_compute_internal_node(unsigned char *		dest,
						     merkle_index_t		node_num,
						     const unsigned char *	seed,
						     param_set_t		lm_type,
						     param_set_t		lm_ots_type,
						     unsigned			h,
						     unsigned			leaf_level,
						     const unsigned char *	I)
{
	int me, proc_num;

	MPI_Comm_rank(MPI_COMM_WORLD, &me);             //进程号
	MPI_Comm_size(MPI_COMM_WORLD, &proc_num);       //进程数
	// struct timespec start, stop;

	// clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start);

	unsigned hash_size = hss_hash_length(h);

	/* We're store intermediate nodes here */
	unsigned char stack[MAX_HASH * MAX_MERKLE_HEIGHT];

	merkle_index_t tree_size = (merkle_index_t)1 << leaf_level;

	merkle_index_t r = node_num;
	unsigned levels_to_bottom = 0;
	//每个进程计算的树高和叶节点数
	unsigned new_levels_to_bottom = 0;
	int leaf_node = 0;

	if (r == 0) return hss_error_internal; /* So no to infinite loops */
	while (r < tree_size) {
		r <<= 1;
		levels_to_bottom++;
	}
	// if(me == 0) printf("me = %d r = %u levels_to_bottom = %u\n",me, r, levels_to_bottom );
	merkle_index_t q = r - tree_size;
	merkle_index_t i;
	unsigned ots_len = lm_ots_get_public_key_len(lm_ots_type);
	unsigned char pub_key[LEAF_MAX_LEN];

	memcpy(pub_key + LEAF_I, I, I_LEN);
	SET_D(pub_key + LEAF_D, D_LEAF);

	struct seed_derive derive;

	if (!hss_seed_derive_init(&derive, lm_type, lm_ots_type,
				  I, seed))
		return hss_error_bad_param_set;

	/*MPI*/

	unsigned char sum_current_buf[proc_num][32];
	int max_proc_num = 0;

	leaf_node = pow(2, levels_to_bottom) / proc_num;
	max_proc_num = leaf_node > 0 ? proc_num : pow(2, levels_to_bottom);
	if (leaf_node > 0) {
		new_levels_to_bottom = log(leaf_node) / log(2);
		r = r + me * leaf_node;
		q = q + me * leaf_node;
	} else if (levels_to_bottom != 0) {
		new_levels_to_bottom = 0;
		r = r + me * 1;
		q = q + me * 1;
	}
	if (me < max_proc_num || levels_to_bottom == 0) {
		for (i = 0;; i++, r++, q++) {
			/* Generate the next OTS public key */

			hss_seed_derive_set_q(&derive, q);
#if (defined(OTS_alltime))
			double start, end;
			start = MPI_Wtime();
#endif
			if (!lm_ots_generate_public_key(lm_ots_type, I,
							q, &derive, pub_key + LEAF_PK, ots_len))
				return hss_error_bad_param_set;
#if (defined(OTS_alltime))
			end = MPI_Wtime();
			sign_kg_count++;
			sign_kg_OTS += (end - start);
#endif
			unsigned char *current_buf;
			int stack_offset = trailing_1_bits(i);
			if (stack_offset == new_levels_to_bottom)
				current_buf = dest;
			else
				current_buf = &stack[stack_offset * hash_size];

			put_bigendian(pub_key + LEAF_R, r, 4);
			union hash_context ctx;
			hss_hash_ctx(current_buf, h, &ctx, pub_key, LEAF_LEN(hash_size));

			unsigned sp;
			for (sp = 1; sp <= stack_offset; sp++) {
				// printf("stack = %02x, current_buf = %02x\n", stack[(sp - 1) * hash_size], current_buf[0]);
				// printf("me = %d h = %u I = %02x hash_size = %u r= %u\n", me, h, I[0], hash_size, r);
				hss_combine_internal_nodes(current_buf,
							   &stack[(sp - 1) * hash_size], current_buf,
							   h, I, hash_size,
							   r >> sp);
				// printf("****r = %u  r >> sp = %u\n",r, r >> sp);
				// printf("****me = %d r >> sp = %u current_buf = %02x\n", me, r >> sp, current_buf[0]);
			}

			if (stack_offset == new_levels_to_bottom)
				break;
		}
	}

	if (levels_to_bottom != 0) {
		// 划分通信域
		// int local_rank[max_proc_num];
		// MPI_Group World_Group, local_group;
		// MPI_Comm local_comm;
		// MPI_Comm_group(MPI_COMM_WORLD, &World_Group);
		//
		// for (int i = 0; i < max_proc_num; i++)
		// 	local_rank[i] = i;
		//
		// MPI_Group_incl(World_Group, max_proc_num, local_rank, &local_group);
		// MPI_Comm_create(MPI_COMM_WORLD, local_group, &local_comm);
		//
		// if (me < max_proc_num)
		// 	MPI_Gather(dest, 32, MPI_CHAR, sum_current_buf[0], 32, MPI_CHAR, 0, local_comm);

		if (me < max_proc_num && me != 0)
			MPI_Send(dest, 32, MPI_CHAR, 0, 0, MPI_COMM_WORLD);

		// if (me == 0) {
		// 	for (size_t i = 1; i < max_proc_num; i++)
		// 		MPI_Recv(sum_current_buf[i], 32, MPI_CHAR, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		// }
		if (me == 0) {
			// MPI_Request request[proc_num - 1];
			// for (size_t i = 1; i < proc_num; i++) {
			// 	MPI_Irecv(sum_current_buf[i], 32, MPI_CHAR, i, 0, MPI_COMM_WORLD, &request[i - 1]);
			// }
			// MPI_Waitall(proc_num - 1, request, MPI_STATUS_IGNORE);
			memcpy(sum_current_buf[0], dest, 32);
			memcpy(dest, "0", 32);

			unsigned levels_to_bottom_0 = levels_to_bottom - new_levels_to_bottom;
			r = r >> new_levels_to_bottom;
			for (i = 0;; i++, r++) {
				if (i != 0)
					MPI_Recv(sum_current_buf[i], 32, MPI_CHAR, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				/* Generate the next OTS public key */

				unsigned char *new_current_buf;
				int stack_offset = trailing_1_bits(i);
				if (stack_offset == levels_to_bottom_0)
					new_current_buf = dest;
				else
					new_current_buf = &stack[stack_offset * hash_size];

				memcpy(new_current_buf, sum_current_buf[i], 32);

				unsigned sp;
				for (sp = 1; sp <= stack_offset; sp++) {
					hss_combine_internal_nodes(new_current_buf,
								   &stack[(sp - 1) * hash_size], new_current_buf,
								   h, I, hash_size,
								   r >> sp);
				}

				if (stack_offset == levels_to_bottom_0)
					break;
			}
		}
		MPI_Bcast(dest, 32, MPI_CHAR, 0, MPI_COMM_WORLD);
	}

	// printf("IIIIIIIIIIIIII\n" );
	// int aaa = 1;
	//
	// while (1)
	// 	aaa++;

	hss_seed_derive_done(&derive);

	// clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &stop);
	// result += (stop.tv_sec - start.tv_sec) * 1e6 + (stop.tv_nsec - start.tv_nsec) / 1e3;
	// if (me == 0) printf("generate the nodee %.2lf msec\n", result / 1e3);

	return hss_error_none;
}
#else /*密钥对生成和签名生成没有并行树*/
enum hss_error_code original_hss_compute_internal_node(unsigned char *		dest,
						       merkle_index_t		node_num,
						       const unsigned char *	seed,
						       param_set_t		lm_type,
						       param_set_t		lm_ots_type,
						       unsigned			h,
						       unsigned			leaf_level,
						       const unsigned char *	I)
{
	unsigned hash_size = hss_hash_length(h);

	/* We're store intermediate nodes here */
	unsigned char stack[MAX_HASH * MAX_MERKLE_HEIGHT];

	merkle_index_t tree_size = (merkle_index_t)1 << leaf_level;
	merkle_index_t r = node_num;
	unsigned levels_to_bottom = 0;

	if (r == 0) return hss_error_internal; /* So no to infinite loops */
	while (r < tree_size) {
		r <<= 1;
		levels_to_bottom++;
	}
	merkle_index_t q = r - tree_size;

	merkle_index_t i;
	unsigned ots_len = lm_ots_get_public_key_len(lm_ots_type);
	unsigned char pub_key[LEAF_MAX_LEN];

	memcpy(pub_key + LEAF_I, I, I_LEN);
	SET_D(pub_key + LEAF_D, D_LEAF);

	struct seed_derive derive;

	if (!hss_seed_derive_init(&derive, lm_type, lm_ots_type,
				  I, seed))
		return hss_error_bad_param_set;

	for (i = 0;; i++, r++, q++) {
		/* Generate the next OTS public key */
		hss_seed_derive_set_q(&derive, q);
#if (defined(MPI_Time)) || (defined(OTS_alltime))
		double start, end;
		start = MPI_Wtime();
#endif
		if (!lm_ots_generate_public_key(lm_ots_type, I,
						q, &derive, pub_key + LEAF_PK, ots_len))
			return hss_error_bad_param_set;         /* The only reason the above */
		                                                /* could fail */

#if (defined(MPI_Time)) || (defined(OTS_alltime))
		end = MPI_Wtime();
		key_count++;
		result_OTS += (end - start);
#endif
		/*
		 * For the subtree which this leaf node forms the final piece, put the
		 * destination to where we'll want it, either on the stack, or if this
		 * is the final piece, to where the caller specified
		 */
		unsigned char *current_buf;
		int stack_offset = trailing_1_bits(i);
		if (stack_offset == levels_to_bottom)
			current_buf = dest;
		else
			current_buf = &stack[stack_offset * hash_size];

		/* Hash it to form the leaf node */
		put_bigendian(pub_key + LEAF_R, r, 4);
		union hash_context ctx;
		hss_hash_ctx(current_buf, h, &ctx, pub_key, LEAF_LEN(hash_size));

		/* Work up the stack, combining right nodes with the left nodes */
		/* that we've already computed */
		unsigned sp;
		for (sp = 1; sp <= stack_offset; sp++) {
			hss_combine_internal_nodes(current_buf,
						   &stack[(sp - 1) * hash_size], current_buf,
						   h, I, hash_size,
						   r >> sp);
		}

		/* We're not at a left branch, or at the target node */

		/* Because we've set current_buf to point to where we want to place */
		/* the result of this loop, we don't need to memcpy it */

		/* Check if this was the last leaf (and so we've just computed the */
		/* target node) */
		if (stack_offset == levels_to_bottom) {
			/* We're at the target node; the node we were asked to compute */
			/* We've already placed the value into dest, so we're all done */
			break;
		}
	}

	hss_seed_derive_done(&derive);

	return hss_error_none;
}

static enum hss_error_code hss_compute_internal_node(unsigned char *		dest,
						     merkle_index_t		node_num,
						     const unsigned char *	seed,
						     param_set_t		lm_type,
						     param_set_t		lm_ots_type,
						     unsigned			h,
						     unsigned			leaf_level,
						     const unsigned char *	I)
{
	unsigned hash_size = hss_hash_length(h);

	/* We're store intermediate nodes here */
	unsigned char stack[MAX_HASH * MAX_MERKLE_HEIGHT];

	merkle_index_t tree_size = (merkle_index_t)1 << leaf_level;
	merkle_index_t r = node_num;
	unsigned levels_to_bottom = 0;

	if (r == 0) return hss_error_internal; /* So no to infinite loops */
	while (r < tree_size) {
		r <<= 1;
		levels_to_bottom++;
	}
	merkle_index_t q = r - tree_size;

	merkle_index_t i;
	unsigned ots_len = lm_ots_get_public_key_len(lm_ots_type);
	unsigned char pub_key[LEAF_MAX_LEN];

	memcpy(pub_key + LEAF_I, I, I_LEN);
	SET_D(pub_key + LEAF_D, D_LEAF);

	struct seed_derive derive;

	if (!hss_seed_derive_init(&derive, lm_type, lm_ots_type,
				  I, seed))
		return hss_error_bad_param_set;

	for (i = 0;; i++, r++, q++) {
		/* Generate the next OTS public key */
		hss_seed_derive_set_q(&derive, q);
#if (defined(OTS_alltime))
		double start, end;
		start = MPI_Wtime();
#endif
		if (!lm_ots_generate_public_key(lm_ots_type, I,
						q, &derive, pub_key + LEAF_PK, ots_len))
			return hss_error_bad_param_set;         /* The only reason the above */
		                                                /* could fail */
#if (defined(OTS_alltime))
		end = MPI_Wtime();
		sign_kg_count++;
		sign_kg_OTS += (end - start);
#endif


		/*
		 * For the subtree which this leaf node forms the final piece, put the
		 * destination to where we'll want it, either on the stack, or if this
		 * is the final piece, to where the caller specified
		 */
		unsigned char *current_buf;
		int stack_offset = trailing_1_bits(i);
		if (stack_offset == levels_to_bottom)
			current_buf = dest;
		else
			current_buf = &stack[stack_offset * hash_size];

		/* Hash it to form the leaf node */
		put_bigendian(pub_key + LEAF_R, r, 4);
		union hash_context ctx;
		hss_hash_ctx(current_buf, h, &ctx, pub_key, LEAF_LEN(hash_size));

		/* Work up the stack, combining right nodes with the left nodes */
		/* that we've already computed */
		unsigned sp;
		for (sp = 1; sp <= stack_offset; sp++) {
			hss_combine_internal_nodes(current_buf,
						   &stack[(sp - 1) * hash_size], current_buf,
						   h, I, hash_size,
						   r >> sp);
		}

		/* We're not at a left branch, or at the target node */

		/* Because we've set current_buf to point to where we want to place */
		/* the result of this loop, we don't need to memcpy it */

		/* Check if this was the last leaf (and so we've just computed the */
		/* target node) */
		if (stack_offset == levels_to_bottom) {
			/* We're at the target node; the node we were asked to compute */
			/* We've already placed the value into dest, so we're all done */
			break;
		}
	}

	hss_seed_derive_done(&derive);

	return hss_error_none;
}
#endif





/*
 * Combine adjacent left and right nodes within the Merkle tree
 * together
 */
void hss_combine_internal_nodes(unsigned char *dest,
				const unsigned char *left_node, const unsigned char *right_node,
				int h, const unsigned char *I, unsigned hash_size,
				merkle_index_t node_num)
{
	unsigned char hash_val[INTR_MAX_LEN];

	memcpy(hash_val + INTR_I, I, I_LEN);
	put_bigendian(hash_val + INTR_R, node_num, 4);
	SET_D(hash_val + INTR_D, D_INTR);

	memcpy(hash_val + INTR_PK, left_node, hash_size);
	memcpy(hash_val + INTR_PK + hash_size, right_node, hash_size);

	// printf("le = %02x %02x, ri = %02x, %02x\n", left_node[0], left_node[1], right_node[0], right_node[1]);
	// for (size_t i = 0; i < INTR_MAX_LEN; i++) {
	// 	printf("%02x", hash_val[i]);
	// 	if (i % 8 == 0) printf("\n");
	// }
	// printf("\n");
	union hash_context ctx;

	hss_hash_ctx(dest, h, &ctx, hash_val, INTR_LEN(hash_size));
}

/*
 * This computes an array of intermediate Merkle nodes given by data
 * This may be run in a worker (non-main) thread
 */
//keygen version
void original_hss_gen_intermediate_tree(const void *data)
{
	const struct intermed_tree_detail *d = data;
	unsigned hash_len = hss_hash_length(d->h);
	unsigned i;

	// #pragma omp parallel for
	// #pragma omp parallel
	// {
	for (i = 0; i < d->node_count; i++) {
		// int my_id = omp_get_thread_num();
		// printf("my_id = %d %d\n", my_id, i);
		unsigned char result[MAX_HASH];
		enum hss_error_code status = original_hss_compute_internal_node(result,
										d->node_num + i,
										d->seed,
										d->lm_type,
										d->lm_ots_type,
										d->h,
										d->tree_height,
										d->I);
		/* Report the results */
		//if (status == hss_error_none) {
		/* Copy out the resulting hash */
		memcpy(d->dest + i * hash_len, result, hash_len);
		//} else {
		//	/* Something went wrong; report the bad news */
		//	*d->got_error = status;
		//	return;
		//}
	}
	// }
}

//sign version
void hss_gen_intermediate_tree(const void *data)
{
	int me, proc_num;

	MPI_Comm_rank(MPI_COMM_WORLD, &me);             //进程号
	MPI_Comm_size(MPI_COMM_WORLD, &proc_num);       //进程数

	const struct intermed_tree_detail *d = data;
	unsigned hash_len = hss_hash_length(d->h);
	unsigned i;

	for (i = 0; i < d->node_count; i++) {
		unsigned char result[MAX_HASH];
		// printf("******************\n");
		enum hss_error_code status = hss_compute_internal_node(result,
								       d->node_num + i,
								       d->seed,
								       d->lm_type,
								       d->lm_ots_type,
								       d->h,
								       d->tree_height,
								       d->I);
		// printf("******************\n");
		/* Report the results */
		if (status == hss_error_none) {
			/* Copy out the resulting hash */
			memcpy(d->dest + i * hash_len, result, hash_len);
		} else {
			/* Something went wrong; report the bad news */
			*d->got_error = status;
			return;
		}
	}
}
