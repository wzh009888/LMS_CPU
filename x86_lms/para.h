// #define PARAM_SHA256 "10/4"
// #define PARAM_SHA256 "15/4"
// #define PARAM_SHA256 "20/4"
// #define PARAM_SHA256 "25/4"
#define PARAM_SHA256 "10/4,10/4"
// #define PARAM_SHA256 "10/4,5/4,5/4"
// #define PARAM_SHA256 "15/4,5/4"
// #define PARAM_SHA256 "10/4,10/4,10/4,10/4"
// #define PARAM_SHA256 "15/4,10/4,10/4,5/4"
// #define PARAM_SHA256 "20/4,10/4,5/4,5/4"
// #define PARAM_SHA256 "20/4,10/4,10/4"

// #define PARAM_SHA256 "5/4,5/4"
// #define PARAM_SHA256 "5/4,5/4,5/4,5/4"
// #define PARAM_SHA256 "5/8,5/8,5/8,5/8"
// #define PARAM_SHA256 "10/8,10/8"
// #define PARAM_SHA256 "10/8,5/8,5/8"
#define ITER_NUM 100
#define MF 2700000000

#define HEIGHT_SUBTREE 2

// #define DP // undefine this to use AP

#ifdef DP

#define USE_OPENSSL 0   /* We use the OpenSSL implementation for SHA-256 */
                         /* (which is quite a bit faster than our portable */
                         /* C version) */
#define vectorization

#else

#define USE_OPENSSL 0

#define parallel_tree   /*只并行分支节点，叶节点不优化*/

#define vectorization   /*在parallel_tree的基础上，使用向量化优化叶节点*/

// #define three_parallel  /*在parallel_tree的基础上，使用多进程+向量化优化叶节点*/
// #define MPI_OTS    /*使用多进程优化叶节点*/
#define vector_512 //默认256位

/*用来测密钥对生成最短时间*/
// #define Gather_mintime//用于测试密钥对生成通信的最短时间，正常使用#define Gather来保证正确性
// #define Igather_mintime
// #define Send_Recv_mintime
// #define Isend_Irecv_mintime

/*用来控制六个通信阶段的最优通信方式，选择全是gather的方式是未通信优化的版本*/
/*01:KG-LMS-Tree  02:lmots-pkgen  03:SIG-LMS-Tree
 * 04:lmots-sign  05:lmots-verify 06：SIG-lmots-pkgen 07：lms_sign*/

#define Gather01
// #define Igather01
// #define Send_Recv01
// #define Isend_Irecv01

/*02:lmots-pkgen 进程少时不调用,进程多时Isend_recv效果最好*/
// #define Gather02
// #define Igather02
// #define Send_Recv02
#define Isend_Irecv02

/*03:SIG-LMS-Tree 进程少或多时send_recv效果最好*/
// #define Gather03
// #define Igather03
#define Send_Recv03
// #define Isend_Irecv03

/*04:lmots-sign 进程少或多时Gather效果最好*/
#define Gather04
// #define Igather04
// #define Send_Recv04
// #define Isend_Irecv04

/*05:lmots-verify 进程少或多时Gather效果最好*/
#define Gather05
// #define Igather05
// #define Send_Recv05
// #define Isend_Irecv05

/*06：SIG-lmots-pkgen*/
// #define Gather06
#define Send_Recv06

/*07：lms_sign 进程少或多时Gather效果最好*/
#define Gather07
// #define Igather07
// #define Send_Recv07
// #define Isend_Irecv07

/*用来控制测试时间输出的*/
#define MPI_Time /*用来测执行一个OTS的操作时间，例如，ots-Keygen，ots-sign，ots-verify。*/
// #define OTS_alltime /*用来测三个子算法中OTS的总时间*/
// #define Communication_time /*用来测通信时间*/

#endif
