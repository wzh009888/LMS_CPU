#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "common_defs.h"
#include "hss.h"
#include "hss_internal.h"
#include "hss_aux.h"
#include "endian.h"
#include "hash.h"
#include "hss_thread.h"
#include "lm_common.h"
#include "lm_ots_common.h"
#include "mpi.h"
#include "hss_derive.h"

/* Count the number of 1 bits at the end (lsbits) of the integer */
/* Do it in the obvious way; straightline code may be faster (no */
/* unpredictable jumps, which are costly), but that would be less scrutable */
static int trailing_1_bits(merkle_index_t n)
{
	int i;

	for (i = 0; n& 1; n >>= 1, i++)
		;
	return i;
}

/*
 * This creates a private key (and the correspond public key, and optionally
 * the aux data for that key)
 * Parameters:
 * generate_random - the function to be called to generate randomness.  This
 *       is assumed to be a pointer to a cryptographically secure rng,
 *       otherwise all security is lost.  This function is expected to fill
 *       output with 'length' uniformly distributed bits, and return 1 on
 *       success, 0 if something went wrong
 * levels - the number of levels for the key pair (2-8)
 * lm_type - an array of the LM registry entries for the various levels;
 *      entry 0 is the topmost
 * lm_ots_type - an array of the LM-OTS registry entries for the various
 *      levels; again, entry 0 is the topmost
 * update_private_key, context - the function that is called when the
 *      private key is generated; it is expected to store it to secure NVRAM
 *      If this is NULL, then the context pointer is reinterpretted to mean
 *      where in RAM the private key is expected to be placed
 * public_key - where to store the public key
 * len_public_key - length of the above buffer; see hss_get_public_key_len
 *      if you need a hint.
 * aux_data - where to store the optional aux data.  This is not required, but
 *      if provided, can be used to speed up the hss_generate_working_key
 *      process;
 * len_aux_data - the length of the above buffer.  This is not fixed length;
 *      the function will run different time/memory trade-offs based on the
 *      length provided
 *
 * This returns true on success, false on failure
 */
int global_level_nodes;
MPI_Comm local_comm;
extern double result_OTS;
extern double key_count;
double sum_time = 0;
/*MPI version*/
bool hss_generate_private_key(
	bool ( *generate_random )(void *output, size_t length),
	unsigned levels,
	const param_set_t *lm_type,
	const param_set_t *lm_ots_type,
	bool (*update_private_key)(unsigned char *private_key,
				   size_t len_private_key, void *context),
	void *context,
	unsigned char *public_key, size_t len_public_key,
	unsigned char *aux_data, size_t len_aux_data,
	struct hss_extra_info *info)
{
	int me, proc_num;

	MPI_Comm_rank(MPI_COMM_WORLD, &me);             //进程号
	MPI_Comm_size(MPI_COMM_WORLD, &proc_num);       //进程数


	struct hss_extra_info info_temp = { 0 };

	if (!info) info = &info_temp;

	if (!generate_random) {
		/* We *really* need random numbers */
		info->error_code = hss_error_no_randomness;
		return false;
	}
	if (levels < MIN_HSS_LEVELS || levels > MAX_HSS_LEVELS) {
		/* parameter out of range */
		info->error_code = hss_error_bad_param_set;
		return false;
	}

	unsigned h0;            /* The height of the root tree */
	unsigned h;             /* The hash function used */
	unsigned size_hash;     /* The size of each hash that would appear in the */

	/* aux data */
	if (!lm_look_up_parameter_set(lm_type[0], &h, &size_hash, &h0)) {
		info->error_code = hss_error_bad_param_set;
		return false;
	}

	/* Check the public_key_len */
	if (4 + 4 + 4 + I_LEN + size_hash > len_public_key) {
		info->error_code = hss_error_buffer_overflow;
		/* public key won't fit in the buffer we're given */
		return false;
	}

	/* If you provide an aux_data buffer, we have to write something */
	/* into it (at least, enough to mark it as 'we're not really using */
	/* aux data) */
	if (aux_data && len_aux_data == 0) {
		/* not enough aux data buffer to mark it as 'not really used' */
		info->error_code = hss_error_bad_aux;
		return false;
	}

	unsigned len_ots_pub = lm_ots_get_public_key_len(lm_ots_type[0]);

	if (len_ots_pub == 0) {
		info->error_code = hss_error_bad_param_set;
		return false;
	}

	unsigned char private_key[PRIVATE_KEY_LEN];

	/* First step: format the private key */
	put_bigendian(private_key + PRIVATE_KEY_INDEX, 0,
		      PRIVATE_KEY_INDEX_LEN);
	if (!hss_compress_param_set(private_key + PRIVATE_KEY_PARAM_SET,
				    levels, lm_type, lm_ots_type,
				    PRIVATE_KEY_PARAM_SET_LEN)) {
		info->error_code = hss_error_bad_param_set;
		return false;
	}
	if (!(*generate_random)(private_key + PRIVATE_KEY_SEED,
				PRIVATE_KEY_SEED_LEN)) {
		info->error_code = hss_error_bad_randomness;
		return false;
	}

	/* Now make sure that the private key is written to NVRAM */
	if (update_private_key) {
		if (!(*update_private_key)(private_key, PRIVATE_KEY_LEN, context)) {
			/* initial write of private key didn't take */
			info->error_code = hss_error_private_key_write_failed;
			hss_zeroize(private_key, sizeof private_key);
			return false;
		}
	} else {
		if (context == 0) {
			/* We weren't given anywhere to place the private key */
			info->error_code = hss_error_no_private_buffer;
			hss_zeroize(private_key, sizeof private_key);
			return false;
		}
		memcpy(context, private_key, PRIVATE_KEY_LEN);
	}

	/* Figure out what would be the best trade-off for the aux level */
	struct expanded_aux_data *expanded_aux_data = 0, aux_data_storage;

	if (aux_data != NULL) {
		aux_level_t aux_level = hss_optimal_aux_level(len_aux_data, lm_type,
							      lm_ots_type, NULL);
		hss_store_aux_marker(aux_data, aux_level);
		// printf("len_aux_data %u\n", len_aux_data);
		// for (size_t i = 0; i < len_aux_data; i++) {
		//     printf("aux = %02x\n", aux_data[i]);
		// }

		/* Set up the aux data pointers */
		expanded_aux_data = hss_expand_aux_data(aux_data, len_aux_data,
							&aux_data_storage, size_hash, 0);
	}

	unsigned char I[I_LEN];
	unsigned char seed[SEED_LEN];

	/*生成根种子和I值的内部函数(基于私有种子)。
	 * 我们这样做(而不是随机选择seed和I)，因此我们不需要将它存储在我们的私钥中;我们可以重新计算它们*/
	if (!hss_generate_root_seed_I_value(seed, I, private_key + PRIVATE_KEY_SEED)) {
		info->error_code = hss_error_internal;
		hss_zeroize(private_key, sizeof private_key);
		return false;
	}

	/* Now, it's time to generate the public key, which means we need to */
	/* compute the entire top level Merkle tree */

	/* First of all, figure out the appropriate level to compute up to */
	/* in parallel.  We'll do the lower of the bottom-most level that */
	/* appears in the aux data, and 4*log2 of the number of core we have */
	unsigned num_cores = 1;
	unsigned level;
	unsigned char *dest = 0;        /* The area we actually write to */
	void *temp_buffer = 0;          /* The buffer we need to free when done */

	for (level = h0 - 1; level > 2; level--) {
		/* If our bottom-most aux data is at this level, we want it */
		// printf("level = %u\n", level);
		if (expanded_aux_data) {
			// printf("AAAAAAAAAAAAAA\n" );
		}
		if (expanded_aux_data->data[level]) {
			// printf("BBBBBBBBBBBBBB\n" );
		}
		if (expanded_aux_data && expanded_aux_data->data[level]) {
			/* Write directly into the aux area */
			dest = expanded_aux_data->data[level];
			break;
		}

		/* If going to a higher levels would mean that we wouldn't */
		/* effectively use all the cores we have, use this level */
		if ((1 << level) < 4 * num_cores) {
			/* We'll write into a temp area; malloc the space */
			size_t temp_buffer_size = (size_t)size_hash << level;
			temp_buffer = malloc(temp_buffer_size);
			if (!temp_buffer)
				/* Couldn't malloc it; try again with s smaller buffer */
				continue;
			/* Use this buffer */
			dest = temp_buffer;
			break;
		}
	}

	// if (me == 0) printf("level = %u\n", level);


	/* Worse comes the worse, if we can't malloc anything, use a */
	/* small backup buffer */
	unsigned char worse_case_buffer[4 * MAX_HASH];

	if (!dest)
		dest = worse_case_buffer;
	/* level == 2 if we reach here, so the buffer is big enough */

	/*
	 * Now, issue all the work items to generate the intermediate hashes
	 * These intermediate passes are potentially computed in parallel;
	 * allowing that is why we use this funky thread_collection and details
	 * structure
	 */
	//struct thread_collection *col = hss_thread_init(info->num_threads);

	struct intermed_tree_detail details;

	/* Set the values in the details structure that are constant */
	details.seed = seed;
	details.lm_type = lm_type[0];
	details.lm_ots_type = lm_ots_type[0];
	details.h = h;
	details.tree_height = h0;
	details.I = I;
	enum hss_error_code got_error = hss_error_none; /* This flag is set */

	/* on an error */
	details.got_error = &got_error;

	merkle_index_t j;
	/* # of nodes at this level */
	merkle_index_t level_nodes = (merkle_index_t)1 << level;

	/* the index of the node we're generating right now */
	merkle_index_t node_num = level_nodes;

	//printf("node_num = %ld\n", node_num);
	/*
	 * We'd prefer not to issue a separate work item for every node; we
	 * might be doing millions of node (if we have a large aux data space)
	 * and we end up malloc'ing a large structure for every work order.
	 * So, if we do have a large number of requires, aggregate them
	 */
	merkle_index_t increment = level_nodes / (10 * num_cores);

	//printf("increment = %ld\n", increment);

#define MAX_INCREMENT 20000
	if (increment > MAX_INCREMENT) increment = MAX_INCREMENT;
	if (increment == 0) increment = 1;

#if (defined(parallel_tree))
	global_level_nodes = level_nodes;
#if (defined(Gather02)) || (defined(Igather02))
	// printf("global_level_nodes = %d\n", global_level_nodes);
	if (proc_num > level_nodes) {
		int new_size;

		new_size = proc_num / global_level_nodes;

		int local_rank[new_size];           //子通信域中的进程数
		MPI_Group World_Group, local_group;

		// MPI_Comm local_comm;

		MPI_Comm_group(MPI_COMM_WORLD, &World_Group);

		//8是level_nodes
		for (int i = me % level_nodes, r = 0; i < proc_num; r++, i += level_nodes) {
			local_rank[r] = i;         //存全局进程号
		}
		MPI_Group_incl(World_Group, new_size, local_rank, &local_group);
		MPI_Comm_create(MPI_COMM_WORLD, local_group, &local_comm);
	}
#endif

	int id;
	if (proc_num > level_nodes)
		id = me % level_nodes;
	else
		id = me;
	size_t num; //当前进程计算的叶节点数
	//level_nodes = (merkle_index_t)1 << 7;
	int num1 = level_nodes / proc_num;
	int num2 = level_nodes % proc_num;
	int local = num1 + ((id < num2)? 1: 0);
	int offset = id * num1 + ((id < num2)? id: num2);
	num = local;

	//MPI_Barrier(MPI_COMM_WORLD); // 同步

	int dest_offset = offset * size_hash;
	details.dest = dest + dest_offset;
	details.node_num = node_num + offset;
	details.node_count = num;

	original_hss_gen_intermediate_tree(&details);

	// printf("me = %d dest = %02x\n", me, dest[0]);
#ifdef OTS_alltime
	if (me == 0)
		printf("result_OTS = %.4lf\n", result_OTS * 1e3);

#endif
#ifdef MPI_Time
	if (me == 0)
		printf("key_count = %lf result_OTS = %.4lf msec %.4lf msec\n", key_count, result_OTS * 1e3, result_OTS * 1e3 / key_count);

#endif
	// MPI_Barrier(MPI_COMM_WORLD); // 同步

#if (defined(Gather_mintime))
/*阻塞集合通信*/
	int r_counts[proc_num];
	int r_dis[proc_num]; //存起始位置

	for (int i = 0; i < num2; i++) {
		r_counts[i] = (num1 + 1) * size_hash;
		r_dis[i] = (i * num1 + i) * size_hash;
	}
	for (int i = num2; i < proc_num; i++) {
		r_counts[i] = num1 * size_hash;
		r_dis[i] = (i * num1 + num2) * size_hash;
	}
	MPI_Comm local_comm002;
//划分通信域，每level_nodes一组
	if (proc_num > level_nodes) {
		int new_size;
		new_size = level_nodes;
		int local_rank002[new_size];           //子通信域中的进程数
		MPI_Group World_Group, local_group;
		MPI_Comm_group(MPI_COMM_WORLD, &World_Group);
		//8是level_nodes
		for (int i = me / level_nodes * level_nodes, r = 0; r < level_nodes; r++, i++) {
			local_rank002[r] = i;         //存全局进程号
		}
		MPI_Group_incl(World_Group, new_size, local_rank002, &local_group);
		MPI_Comm_create(MPI_COMM_WORLD, local_group, &local_comm002);
	}
//end 划分通信域
#if (defined(Communication_time))
	MPI_Barrier(MPI_COMM_WORLD); // 同步
	double t0, t1, t2, t_sum;
	t0 = MPI_Wtime();
#endif

	if (proc_num > level_nodes)
		MPI_Gatherv(&dest[dest_offset], num * size_hash, MPI_CHAR, dest, r_counts, r_dis, MPI_CHAR, 0, local_comm002);
	else
		MPI_Gatherv(&dest[dest_offset], num * size_hash, MPI_CHAR, dest, r_counts, r_dis, MPI_CHAR, 0, MPI_COMM_WORLD);

#if (defined(Communication_time))
	t1 = MPI_Wtime();
	t2 = t1 - t0;
// printf("t2 = %f\n", t2 * 1000);

	MPI_Reduce(&t2, &t_sum, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	sum_time += (t_sum * 1000);
	if (me == 0)
		printf("************sum LMS-Tree of KG Gather_time = %.6lf***************\n", sum_time);
#endif

#elif (defined(Igather_mintime))
	/*非阻塞集合通信*/
	int r_counts[proc_num];
	int r_dis[proc_num];        //存起始位置

	for (int i = 0; i < num2; i++) {
		r_counts[i] = (num1 + 1) * size_hash;
		r_dis[i] = (i * num1 + i) * size_hash;
	}
	for (int i = num2; i < proc_num; i++) {
		r_counts[i] = num1 * size_hash;
		r_dis[i] = (i * num1 + num2) * size_hash;
	}
	MPI_Comm local_comm002;
	//划分通信域，每level_nodes一组
	if (proc_num > level_nodes) {
		int new_size;
		new_size = level_nodes;
		int local_rank002[new_size];                   //子通信域中的进程数
		MPI_Group World_Group, local_group;
		MPI_Comm_group(MPI_COMM_WORLD, &World_Group);
		//8是level_nodes
		for (int i = me / level_nodes * level_nodes, r = 0; r < level_nodes; r++, i++) {
			local_rank002[r] = i;                 //存全局进程号
		}
		MPI_Group_incl(World_Group, new_size, local_rank002, &local_group);
		MPI_Comm_create(MPI_COMM_WORLD, local_group, &local_comm002);
	}
	//end 划分通信域
#if (defined(Communication_time))
	MPI_Barrier(MPI_COMM_WORLD);         // 同步
	double t0, t1, t2, t_sum;
	t0 = MPI_Wtime();
#endif
	MPI_Request re;
	// printf("dest_offset = %d %d\n", dest_offset, me);
	if (proc_num > level_nodes)
		MPI_Igatherv(&dest[dest_offset], num * size_hash, MPI_CHAR, dest, r_counts, r_dis, MPI_CHAR, 0, local_comm002, &re);
	else
		MPI_Igatherv(&dest[dest_offset], num * size_hash, MPI_CHAR, dest, r_counts, r_dis, MPI_CHAR, 0, MPI_COMM_WORLD, &re);
	MPI_Wait(&re, MPI_STATUS_IGNORE);
#if (defined(Communication_time))
	t1 = MPI_Wtime();
	t2 = t1 - t0;
	// printf("t2 = %f\n", t2 * 1000);

	MPI_Reduce(&t2, &t_sum, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	sum_time += (t_sum * 1000);
	if (me == 0)
		printf("************sum LMS-Tree of KG Igather_time = %.6lf***************\n", sum_time);
#endif

#elif (defined(Send_Recv_mintime))
#if (defined(Communication_time))
	MPI_Barrier(MPI_COMM_WORLD);         // 同步
	//阻塞点对点通信所有进程将数据发送给0号进程，0号进程将收集完的数据发送给其他进程
	double t0, t1, t2, t_sum;
	t0 = MPI_Wtime();
#endif
	if (me != 0)
		MPI_Send(&dest[dest_offset], num * size_hash, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
	if (me == 0) {
		if (num2 > 1) {
			for (size_t i = 1; i < num2; i++) {
				unsigned int b = (i * num1 + i) * size_hash;
				MPI_Recv(&dest[b], (num1 + 1) * size_hash, MPI_CHAR, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			}
			for (size_t i = num2; i < proc_num; i++) {
				unsigned int b = (i * num1 + num2) * size_hash;
				MPI_Recv(&dest[b], num1 * size_hash, MPI_CHAR, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			}
		} else {
			for (size_t i = 1; i < proc_num; i++) {
				unsigned int b = (i * num1 + num2) * size_hash;
				MPI_Recv(&dest[b], num1 * size_hash, MPI_CHAR, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			}
		}
	}

#if (defined(Communication_time))
	t1 = MPI_Wtime();
	t2 = t1 - t0;
	// printf("t2 = %f\n", t2 * 1000);

	MPI_Reduce(&t2, &t_sum, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	sum_time += (t_sum * 1000);
	if (me == 0)
		printf("************sum LMS-Tree of KG Send_Recv_time = %.6lf***************\n", sum_time);
#endif

#elif (defined(Isend_Irecv_mintime))
	//非阻塞点对点通信，所有进程将自己的数据发送给其他所有进程，再接收其他所有进程发送给自己的数据
	//发送数据
	MPI_Comm local_comm002;
	//划分通信域，每level_nodes一组
	if (proc_num > level_nodes) {
		int new_size;
		new_size = level_nodes;
		int local_rank002[new_size];                   //子通信域中的进程数
		MPI_Group World_Group, local_group;
		MPI_Comm_group(MPI_COMM_WORLD, &World_Group);

		for (int i = me / level_nodes * level_nodes, r = 0; r < level_nodes; r++, i++) {
			local_rank002[r] = i;                 //存全局进程号
		}
		MPI_Group_incl(World_Group, new_size, local_rank002, &local_group);
		MPI_Comm_create(MPI_COMM_WORLD, local_group, &local_comm002);
	}
	//end 划分通信域
	unsigned int sub_num = 0;
	if (proc_num > level_nodes) sub_num = level_nodes;
	else sub_num = proc_num;
#if (defined(Communication_time))
	MPI_Barrier(MPI_COMM_WORLD);         // 同步
	double t0, t1, t2, t_sum;
	t0 = MPI_Wtime();
#endif
	if (proc_num < level_nodes) {         //即一个进程算多个分支节点的值
		MPI_Request request1[proc_num], request2[sub_num];
		for (size_t i = 0; i < proc_num; i++)
			MPI_Isend(&dest[dest_offset], num * size_hash, MPI_CHAR, i, 0, MPI_COMM_WORLD, &request1[i]);
		for (size_t j = 0; j < proc_num; j++) {
			unsigned int b = (j * num1) * size_hash;
			MPI_Irecv(&dest[b], num1 * size_hash, MPI_CHAR, j, 0, MPI_COMM_WORLD, &request2[j]);
		}
		MPI_Waitall(proc_num, request1, MPI_STATUS_IGNORE);
		MPI_Waitall(sub_num, request2, MPI_STATUS_IGNORE);
	} else {
		MPI_Request request1[level_nodes], request2[level_nodes];
		// if (me < 8) {
		for (size_t i = 0; i < level_nodes; i++)
			MPI_Isend(&dest[dest_offset], num * size_hash, MPI_CHAR, i, 0, local_comm002, &request1[i]);
		for (size_t j = 0; j < level_nodes; j++)
			MPI_Irecv(&dest[j * size_hash], num * size_hash, MPI_CHAR, j, 0, local_comm002, &request2[j]);
		int me1, proc_num1;
		if (me < 16 && me >= 8) {
			MPI_Comm_rank(local_comm002, &me1);                     //进程号
			MPI_Comm_size(local_comm002, &proc_num1);               //进程数
			// printf("!!!! %d %d %d %d %d\n", me1, proc_num1, num, size_hash, dest_offset);
		}
		// printf("me = %d, b = %d %d\n", me, num1, size_hash);
		MPI_Waitall(level_nodes, request1, MPI_STATUS_IGNORE);
		MPI_Waitall(level_nodes, request2, MPI_STATUS_IGNORE);
	}
#if (defined(Communication_time))
	t1 = MPI_Wtime();
	t2 = t1 - t0;
	// printf("t2 = %f\n", t2 * 1000);

	MPI_Reduce(&t2, &t_sum, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	sum_time += (t_sum * 1000);
	if (me == 0)
		printf("************sum LMS-Tree of KG Isend_Irecv_time = %.6lf***************\n", sum_time);
#endif
#endif
	hss_zeroize(seed, sizeof seed);

#else
	for (j = 0; j < level_nodes;) {
		unsigned this_increment;
		if (level_nodes - j < increment)
			this_increment = level_nodes - j;
		else
			this_increment = increment;

		/* Set the particulars of this specific work item */
		details.dest = dest + j * size_hash;
		details.node_num = node_num;
		details.node_count = this_increment;

		/* Issue a separate work request for every node at this level */
		original_hss_gen_intermediate_tree(&details);

		j += this_increment;
		node_num += this_increment;
	}
#ifdef OTS_alltime
	if (me == 0)
		printf("result_OTS = %.4lf\n", result_OTS * 1e3);

#endif
#ifdef MPI_Time
	if (me == 0)
		printf("key_count = %lf Original_result_OTS = %.4lf msec %.4lf msec\n", key_count, result_OTS * 1e3, result_OTS * 1e3 / key_count);

#endif
	hss_zeroize(seed, sizeof seed);
#endif
	/* Check if something went wrong.  It really shouldn't have, however if */
	/* something returns an error code, we really should try to handle it */
	if (got_error != hss_error_none) {
		/* We failed; give up */
		info->error_code = got_error;
		hss_zeroize(private_key, sizeof private_key);
		if (update_private_key)
			(void)(*update_private_key)(private_key, PRIVATE_KEY_LEN, context);
		else
			hss_zeroize(context, PRIVATE_KEY_LEN);
		free(temp_buffer);
		return false;
	}

	/* Now, we complete the rest of the tree.  This is actually fairly fast */
	/* (one hash per node) so we don't bother to parallelize it */

	unsigned char stack[MAX_HASH * (MAX_MERKLE_HEIGHT + 1)];
	unsigned char root_hash[MAX_HASH];

	/* Generate the top levels of the tree, ending with the root node */
	merkle_index_t r, leaf_node;

	// if (me == 0) {
		for (r = level_nodes, leaf_node = 0; leaf_node < level_nodes; r++, leaf_node++) {
			/* Walk up the stack, combining the current node with what's on */
			/* the atack */
			merkle_index_t q = leaf_node;

			/*
			 * For the subtree which this leaf node forms the final piece, put the
			 * destination to where we'll want it, either on the stack, or if this
			 * is the final piece, to where the caller specified
			 */
			unsigned char *current_buf;
			int stack_offset = trailing_1_bits(leaf_node);
			if (stack_offset == level)
				current_buf = root_hash;
			else
				current_buf = &stack[stack_offset * size_hash];
			memcpy(current_buf, dest + leaf_node * size_hash, size_hash);

			unsigned sp;
			unsigned cur_lev = level;
			for (sp = 1;; sp++, cur_lev--, q >>= 1) {
				/* Give the aux data routines a chance to save the */
				/* intermediate value.  Note that we needn't check for the */
				/* bottommost level; if we're saving aux data at that level, */
				/* we've already placed it there */
				if (sp > 1)
					hss_save_aux_data(expanded_aux_data, cur_lev,
							  size_hash, q, current_buf);

				if (sp > stack_offset) break;
				hss_combine_internal_nodes(current_buf,
							   &stack[(sp - 1) * size_hash], current_buf,
							   h, I, size_hash,
							   r >> sp);
			}
		}
		/* The top entry in the stack is the root value (aka the public key) */



		/* Complete the computation of the aux data */
		hss_finalize_aux_data(expanded_aux_data, size_hash, h,
				      private_key + PRIVATE_KEY_SEED);
	// }
	// printf("me = %d\n", me);
	// MPI_Barrier(MPI_COMM_WORLD);


	/* We have the root value; now format the public key */
	put_bigendian(public_key, levels, 4);
	public_key += 4; len_public_key -= 4;
	put_bigendian(public_key, lm_type[0], 4);
	public_key += 4; len_public_key -= 4;
	put_bigendian(public_key, lm_ots_type[0], 4);
	public_key += 4; len_public_key -= 4;
	memcpy(public_key, I, I_LEN);
	public_key += I_LEN; len_public_key -= I_LEN;
	memcpy(public_key, root_hash, size_hash);
	public_key += size_hash; len_public_key -= size_hash;

	/* Hey, what do you know -- it all worked! */
	hss_zeroize(private_key, sizeof private_key);   /* Zeroize local copy of */
	                                                /* the private key */

	free(temp_buffer);
	return true;
}

/*
 * The length of the private key
 */
size_t hss_get_private_key_len(unsigned			levels,
			       const param_set_t *	lm_type,
			       const param_set_t *	lm_ots_type)
{
	/* A private key is a 'public object'?  Yes, in the sense that we */
	/* export it outside this module */
	return PRIVATE_KEY_LEN;
}
