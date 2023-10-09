#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <ctype.h>
#include "mpi.h"
#include "hss.h"
#include "hss_verify_inc.h"
#include "hss_sign_inc.h"
#include "unistd.h"

const char *default_parm_set = PARAM_SHA256;

unsigned char privkey[HSS_MAX_PRIVATE_KEY_LEN];
static size_t got_len;
static unsigned long last_seqno;
static bool got_update;
static bool got_error;
static bool hit_end;
double sum_time6 = 0;
extern double sign_kg_OTS;
extern double sign_kg_count;
extern int main_count;
#define DEFAULT_AUX_DATA 10916   /* Use 10+k of aux data (which works well */
/* with the above default parameter set) */

static const char *seedbits = 0;
static const char *i_value = 0;
typedef unsigned long long u64;

static bool convert_specified_seed_i_value(void *, size_t);

#include "hash.h"
#include "hss_zeroize.h"
// bool rand_1( void *output, size_t len ) {
//     if (seedbits) {
//         /* The seed was specified on the command line */
//         /* Return that exact seed and i */
//         /* This is not something a real application should do */
//         return convert_specified_seed_i_value( output, len );
//     }
//     struct {
//         unsigned char dev_random_output[32];
//         int rand_output[16];
//         /* Potentially more random sources here */
//         unsigned count;
//     } buffer;
//     int i;
//
//     /* Try to grab a sammple of /dev/urandom output */
//     /* We use /dev/urandom because there's no point in blocking; this is a */
//     /* demo program */
//     FILE *f = fopen( "/dev/urandom", "r" );
//     if (f) {
//          (void)fread( buffer.dev_random_output, 1, 32, f );
//          fclose(f);
//     }
//
//     /* Also try to grab some output from rand */
//     /* It's not great, but if the /dev/urandom output fails, at least we */
//     /* have something */
//     /* In a real program, we'd want to fail if we don't have enough */
//     /* entropy, but hey, this is a demo */
//     static int set_seed = 0;
//     if (!set_seed) {
//         srand( time(0) );
//         set_seed = 1;
//     }
//     for (i = 0; i<16; i++) {
//         buffer.rand_output[i] = rand();
//     }
//
//
//     /* If we had more random sources, we'd sample them here */
//
//     unsigned output_buffer[32];
//     for (i=0; len>0; i++) {
//         buffer.count = i;
//
//         /* Ok, hash all our random samples together to generate the random */
//         /* string that was asked for */
//         hss_hash( output_buffer, HASH_SHA256, &buffer, sizeof buffer );
//
//         /* Copy that hash to the output buffer */
//         int this_len = 32;
//         if (this_len > len) this_len = len;
//         memcpy( output, output_buffer, this_len );
//
//         /* Advance pointers */
//         output = (unsigned char *)output + this_len; len -= this_len;
//     }
//
//     /* Clean up after ourselves.  Yes, this is a demo program; doesn't mean */
//     /* we get to be sloppy */
//     hss_zeroize( output_buffer, sizeof output_buffer );
//     hss_zeroize( &buffer, sizeof buffer );
//
//     return true;
// }
static bool rand_1(void *output, size_t len)
{
	got_len = len;
	unsigned char *p = output;

	while (len--) *p++ = 0x01;
	return true;
}


static int get_integer(const char **p)
{
	int n = 0;

	while (isdigit(**p)) {
		n = 10 * n + **p - '0';
		*p += 1;
	}

	return n;
}

/*
 * This parses the parameter set; this is provided so we can try different
 * sets without recompiling the program each time.  This is placed here
 * because it's ugly parsing code that has nothing to do with how to use
 * HSS
 */
static int parse_parm_set(int *levels, param_set_t *lm_array,
			  param_set_t *ots_array, size_t *aux_size,
			  const char *parm_set)
{
	int i;
	size_t aux = DEFAULT_AUX_DATA;

	for (i = 0;; i++) {
		if (i == 8) {
			printf("Error: more than 8 HSS levels specified\n");
			return 0;
		}
		/* Get the number of levels of this tree */
		int h = get_integer(&parm_set);
		param_set_t lm;
		switch (h) {
		case 5:  lm = LMS_SHA256_N32_H5;  break;
		case 10: lm = LMS_SHA256_N32_H10; break;
		case 15: lm = LMS_SHA256_N32_H15; break;
		case 20: lm = LMS_SHA256_N32_H20; break;
		case 25: lm = LMS_SHA256_N32_H25; break;
		case 0: printf("Error: expected height of Merkle tree\n"); return 0;
		default: printf("Error: unsupported Merkle tree height %d\n", h);
			printf("Supported heights = 5, 10, 15, 20, 25\n");
			return 0;
		}
		/* Now see if we can get the Winternitz parameter */
		param_set_t ots = LMOTS_SHA256_N32_W8;
		if (*parm_set == '/') {
			parm_set++;
			int w = get_integer(&parm_set);
			switch (w) {
			case 1: ots = LMOTS_SHA256_N32_W1; break;
			case 2: ots = LMOTS_SHA256_N32_W2; break;
			case 4: ots = LMOTS_SHA256_N32_W4; break;
			case 8: ots = LMOTS_SHA256_N32_W8; break;
			case 0: printf("Error: expected Winternitz parameter\n"); return 0;
			default: printf("Error: unsupported Winternitz parameter %d\n", w);
				printf("Supported parmaeters = 1, 2, 4, 8\n");
				return 0;
			}
		}

		lm_array[i] = lm;
		ots_array[i] = ots;

		if (*parm_set == ':') {
			parm_set++;
			aux = get_integer(&parm_set);
			break;
		}
		if (*parm_set == '\0') break;
		if (*parm_set == ',') {
			parm_set++; continue;
		}
		printf("Error: parse error after tree specification\n"); return 0;
	}

	*levels = i + 1;
	*aux_size = aux;
	return 1;
}

// static bool rand_1(void *output, size_t len) {
//     unsigned char *p = output;
//     while (len--) *p++ = len;
//     return true;
// }


struct priv_key_reader {
	size_t		length;
	unsigned char	priv_key[100];
};

static bool update_private_key(unsigned char *private_key,
			       size_t len_private_key, void *context)
{
	if (len_private_key > HSS_MAX_PRIVATE_KEY_LEN || len_private_key < 8) return false;

	memcpy(privkey, private_key, len_private_key);

	got_update = true;
	hit_end = false;
	got_error = false;

	int i;

	for (i = 0; i < 8; i++)
		if (private_key[i] != 0xff) break;
	if (i == 8) {
		hit_end = true;
		return true;
	}

	/* Our tests never have seqno's larger than 2**32-1 */
	/* If we see any larger claimed, it's an error */
	for (i = 0; i < 4; i++) {
		if (private_key[i] != 0x00) {
			got_error = true;
			return true;
		}
	}

	/* Pull out the sequence number from the private key */
	last_seqno = 0;
	for (i = 4; i < 8; i++)
		last_seqno = 256 * last_seqno + private_key[i];

	return true;
}

static bool read_private_key(unsigned char *private_key,
			     size_t len_private_key, void *context)
{
	if (len_private_key > HSS_MAX_PRIVATE_KEY_LEN) return false;
	memcpy(private_key, privkey, len_private_key);
	return true;
}

u64 cpucycles(void)
{
	u64 result;
	__asm volatile (".byte 15;.byte 49;shlq $32,%%rdx;orq %%rdx,%%rax"
			: "=a" (result) ::  "%rdx");

	return result;
}

/*Create multiple HSS trees*/
static bool data_parallelism_genkey_gensig_valisig(const char *parm_set)
{
	int me, proc_num;

	MPI_Comm_rank(MPI_COMM_WORLD, &me);             //进程号
	MPI_Comm_size(MPI_COMM_WORLD, &proc_num);       //进程数
	unsigned allaux_len[proc_num];
	unsigned allsig_size[proc_num];

	char message[4];
	int levels;
	param_set_t lm_array[MAX_HSS_LEVELS];
	param_set_t ots_array[MAX_HSS_LEVELS];
	size_t aux_size;

	struct timespec start, stop;
	double result;
	u64 t0, t1;
	double mf = MF;
	int inter = ITER_NUM;

	if (!parm_set) parm_set = default_parm_set;
	if (!parse_parm_set(&levels, lm_array, ots_array, &aux_size, parm_set))
		return false;

	int pubkey_size = hss_get_public_key_len(levels, lm_array, ots_array);
	int privkey_size = hss_get_private_key_len(levels, lm_array, ots_array);
	int sig_size = hss_get_signature_len(levels, lm_array, ots_array);

	if (!pubkey_size || pubkey_size > HSS_MAX_PUBLIC_KEY_LEN ||
	    !sig_size ||
	    !privkey_size || privkey_size > HSS_MAX_PRIVATE_KEY_LEN) {
		printf("Internal error: bad parm set\n");
		return false;
	}
	unsigned char pubkey[HSS_MAX_PUBLIC_KEY_LEN];
	unsigned char *sig = malloc(sig_size);

	if (!sig) return false;

	MPI_Gather(&sig_size, 1, MPI_INT, allsig_size, 1, MPI_INT, 0, MPI_COMM_WORLD);


	/* And we'll place the aux data in this array */
	unsigned aux_len;

	if (aux_size > 0) {
		aux_len = hss_get_aux_data_len(aux_size, levels,
					       lm_array, ots_array);
		// printf("aux_len = %d\n", aux_len);
	} else {
		aux_len = 1;
	}
	unsigned char *aux = malloc(aux_len);

	MPI_Gather(&aux_len, 1, MPI_INT, allaux_len, 1, MPI_INT, 0, MPI_COMM_WORLD);


	if (!aux) {
		printf("error mallocing aux; not generating aux\n");
		aux_len = 0;
		aux = 0;
	}

	/*First, generate the key*/
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start);
	t0 = cpucycles();
	struct hss_extra_info info = { 1 };

	for (size_t i = 0; i < inter; i++) {
		/* code */
		if (!hss_generate_private_key(rand_1, levels, lm_array, ots_array,
					      update_private_key, privkey, pubkey, pubkey_size,
					      aux_size > 0 ? aux : 0, aux_len, &info)) {
			printf("Pubkey gen failure\n");
			return false;
		}
	}

	/*Pass the key pair back to process 0*/
	unsigned char allprivkey[proc_num * privkey_size];
	unsigned char allpubkey[proc_num * pubkey_size];
	unsigned char allaux[proc_num * aux_len];
	unsigned char allsig[proc_num * sig_size];
	unsigned char allmess[proc_num * 4];

	if (me == 0)
		for (size_t i = 0; i < proc_num; i += 4)
			memcpy(&allmess[i], "123\n", 4);

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Gather(privkey, privkey_size, MPI_CHAR, allprivkey, privkey_size, MPI_CHAR, 0, MPI_COMM_WORLD);
	MPI_Gather(pubkey, pubkey_size, MPI_CHAR, allpubkey, pubkey_size, MPI_CHAR, 0, MPI_COMM_WORLD);
	MPI_Gather(aux, aux_len, MPI_CHAR, allaux, aux_len, MPI_CHAR, 0, MPI_COMM_WORLD);

	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &stop);
	t1 = cpucycles();
	result = (stop.tv_sec - start.tv_sec) * 1e6 + (stop.tv_nsec - start.tv_nsec) / 1e3;
	if (me == 0) {
		printf("KG %.2lf ms (using time)\n", result / 1e3 / inter);
		printf("KG %.2lf ms (using cycle)\n", (t1 - t0) / mf * 1e3 / inter);
		printf("machine cycle = %.2lld\n", (t1 - t0) / inter);
        printf("KG Throughput = %.3lf\n", proc_num / result * 1e6 * inter);
        printf("-----------------------------\n");
	}

	/*Process 0 sends privkey and messages to other processes*/
	MPI_Barrier(MPI_COMM_WORLD); // 同步
	MPI_Scatter(allprivkey, privkey_size, MPI_CHAR, privkey, privkey_size, MPI_CHAR, 0, MPI_COMM_WORLD);
	MPI_Scatter(allaux, aux_len, MPI_CHAR, aux, aux_len, MPI_CHAR, 0, MPI_COMM_WORLD);
	MPI_Scatter(allaux_len, 1, MPI_INT, &aux_len, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Scatter(allmess, 4, MPI_CHAR, message, 4, MPI_CHAR, 0, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD); // 同步
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start);
	t0 = cpucycles();
	struct hss_working_key *w;

	// message = "Test Messages\n";
	for (size_t i = 0; i < inter; i++) {
		/*Next, generate the signature*/

		w = hss_load_private_key(
			read_private_key, privkey,      /* How to load the */
			/* private key */
			0,                              /* Use minimal memory */
			aux, aux_len,                   /* The auxiliary data */
			&info);                         /* Use the defaults for extra info */

		if (!w) {
			printf("id = %d Error loading working key\n", me);
			return false;
		}


		// struct hss_extra_info info = { 0 };
		bool success = hss_generate_signature(w, update_private_key, privkey,
						      message, sizeof message,
						      sig, sig_size, &info);

		if (!success) {
			printf("Error: Signature failure\n");
			hss_free_working_key(w);
			free(sig);
			return false;
		}
	}


	/*Send the signature to process 0*/
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Gather(privkey, privkey_size, MPI_CHAR, allprivkey, privkey_size, MPI_CHAR, 0, MPI_COMM_WORLD);
	MPI_Gather(sig, sig_size, MPI_CHAR, allsig, sig_size, MPI_CHAR, 0, MPI_COMM_WORLD);

	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &stop);
	t1 = cpucycles();
	result = (stop.tv_sec - start.tv_sec) * 1e6 + (stop.tv_nsec - start.tv_nsec) / 1e3;
	if (me == 0) {
		printf("SIGN %.2lf ms (using time)\n", result / 1e3 / inter);
		printf("SIGN %.2lf ms (using cycle)\n", (t1 - t0) / mf * 1e3 / inter);
		printf("machine cycle =  %.2lld \n", (t1 - t0) / inter);
        printf("SIGN Throughput = %.3lf\n", proc_num / result * 1e6 * inter);
        printf("-----------------------------\n");
	}

	/*Process 0 sends the public key and signature to other processes*/
	MPI_Barrier(MPI_COMM_WORLD);
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start);
	t0 = cpucycles();

	MPI_Scatter(allpubkey, pubkey_size, MPI_CHAR, pubkey, pubkey_size, MPI_CHAR, 0, MPI_COMM_WORLD);
	MPI_Scatter(allsig, sig_size, MPI_CHAR, sig, sig_size, MPI_CHAR, 0, MPI_COMM_WORLD);
	MPI_Scatter(allsig_size, 1, MPI_INT, &sig_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Scatter(allmess, 1, MPI_CHAR, message, 1, MPI_CHAR, 0, MPI_COMM_WORLD);


	/*Third, verify the signature*/
	int ret = 0;

	for (size_t i = 0; i < inter; i++) {
		MPI_Barrier(MPI_COMM_WORLD);
		struct hss_extra_info info = { 0 };
		ret = hss_validate_signature(pubkey, message, sizeof message,
					     sig, sig_size, &info);
		ret += ret;
	}

	if (me == 0) {
		if (ret != 0) {
			printf("VERIFY success\n");
		} else {
			printf("    *** verification failed when it should have passed\n");
			hss_free_working_key(w);
			free(sig);
			return false;
		}
	}

	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &stop);
	t1 = cpucycles();
	result = (stop.tv_sec - start.tv_sec) * 1e6 + (stop.tv_nsec - start.tv_nsec) / 1e3;
	if (me == 0) {
		printf("VERIFY %.2lf ms (using time)\n", result / 1e3 / inter);
		printf("VERIFY %.2lf ms (using cycle)\n", (t1 - t0) / mf * 1e3 / inter);
		printf("machine cycle = %.2lld\n", (t1 - t0) / inter);
        printf("VERIFY Throughput %.3lf tasks/s\n", proc_num / result * 1e6 * inter);
        printf("-----------------------------\n");
	}
	hss_free_working_key(w);

	return true;
}

int main(int argc, char **argv)
{
	MPI_Init(NULL, NULL);
	int me, proc_num;

	MPI_Comm_rank(MPI_COMM_WORLD, &me);             //进程号
	MPI_Comm_size(MPI_COMM_WORLD, &proc_num);       //进程数
	/* Parse the parameter set */
	const char *parm_set = 0;
    if(me == 0) printf("******* data parallelism start ********\n");
    if(me == 0) printf("-----------------------------\n");

    if(me == 0) {
        if(USE_OPENSSL)
            printf("USING OPENSSL\n");
        else
            printf("NOT USING OPENSSL\n");
    }

	if(!data_parallelism_genkey_gensig_valisig(parm_set)){
	    printf("*******dp error********\n");
	}
    MPI_Barrier(MPI_COMM_WORLD);
    if(me == 0) printf("*******dp success********\n");

	MPI_Finalize();

	return 0;
}
