#ifndef _PARAMS_H_
#define _PARAMS_H_

#define PITTPACKACC _OPENACC
#define EXACT 0
#define DEBUG2 0 /*!< turns on debug for poissonGPU */
#define DEBUG_COMM 0
#define DEBUG1 0 /*!<mainly turn off IO for debugging */
#define DEBUG0 0 /*!<mainly turn off IO for debugging  for solverCPU */
#define DEBUG 0  /*!<mainly turn off IO for debugging  for solverCPU */
#define COMM_PATTERN                                                                                                                       \
    2 /*!< two type of comm_patterns are available, (0) Pairwise exchange (1) Neighorhood collective, don not try comm_pattern with CPU    \
         yet */
#define SHORT_                                                                                                                             \
    0 /*!< controls the type for \f$iax, iay, jax, jay\f$ arrays, setting to 1 selects \f$short int\f$  where settingt it to zero will set \
         the type as \f$int\f$ */
#define ZSIZE                                                                                                                              \
    64 /*! shoud be nzChubk*nChunk, not required if the first 4 letters specifying the boundary are not 'P', also not required in new      \
          version of solve  */
#define MULTIGRIDON 0
#define NXCHUNK1 256
#define NYCHUNK1 256
#define TMPSIZE 256

#define COMM_ON 1

#define DISABLE_CHECKS 1

#define MIN( x, y ) ( ( ( x ) < ( y ) ) ? ( x ) : ( y ) )

// const int CHUNK=5;
// Please note that if your ny*nz*sqrt(COM_SIZE) < 65535
// you can use unsigned short instead of int and save memory by the factor of 2
// considering block sizes as multiplis of 16, you should be good to go with short int
// if not go to the file definitions.h and change the typedef there

enum PittPackParams /*!<Parameters to set before compiling  */
{
    I_O                         = 1, /*!< set o zero to disable IO */
    SHIFT                       = 1, /*!< Set to 0 for no shifts and 1 to shift*/
    This                        = 0, /*!<Junk for debugging, will be removed later */
    GPUAW                       = 1, /*! chunkwise send/recieve for ZX rotation */
    GPUAW2                      = 1, /*! chunkwise send/recieve for XY rotation */
    MPI_ERROR_DISABLE           = 0, /*! if set to 1 it will reset MPI_ERROR_FATAl to MPI_ERROR_RETURN */
    PITT_ABORT                  = 0,
    INCLUDE_ERROE_CAL_IN_TIMING = 0, /*! uses an expensive allreduce function misleading to be included in profiling, takes 8 % of 512M
                                        mesh, suggest truning it off for profiling*/
    SOLUTIONMETHOD = 0,              /*!< (0) solves with Thomas (1) Uses PCR (2) CR-P (4) Multigrid and (5) cuSPARCE (CR and PCR )
                                          the last two are disabled disabled to avoid unnecesary memory usage,
                                           need to change the MULTIGRIDON to 1 to enable allocation for multigrid and cuSPARSE */
    INITANALYTIC = 0, /*! This enables initializeTrigonometric function used for debugging and verification of oder of accuracy*/
    PIVOT        = 1, /*!< This is only used for cuSPARSE, for diagonally dominant matrix pivoting is not required */
    INNERITER    = 10,
    OUTERITER    = 20,
    TRI_NUM_GANG
    = 10000, /*!< Overrides the default number of gpu blocks (gangs), setting this to a large value will replace this with nyChunk */
    VECLENGTH     = 256, /*!< sets the number of gpu-threads */
    MAXLEVEL      = 20,
    COEFF0        = 2, /*!< (COEFF0*pi) is the frequency of the exact solution in the x-direction  */
    COEFF1        = 2, /*!< (COEFF0*pi) is the frequency of the exact solution in the y-direction  */
    COEFF2        = 2, /*!< (COEFF0*pi) is the frequency of the exact solution in the z-direction  */
    CUFFT_STREAMS = 0, /*!< uses different stream for fft transform  */
    MONITOR_MEM   = 0, /*!< set this to 1 to monitor memory usage at each level of solution  */
    PROFILE_COMM  = 0, /*!< Turn it on and it will report the amount of time spent for communication */
    NITER         = 1, /*!< Number of iteration for Poisson Solver*/

};

typedef enum PittPackErrorCodes
{
    SUCCESS                        = 0, /*!< Succesful execution */
    GPU_INIT_FAILURE               = 1, /*!< GPU initialization and binding failed*/
    SOLVER_ARRAY_INCONSISTENT      = 2,
    MPI_INIT_FAIL                  = 3, /*!< MPI initialization failed */
    MPI_DUP_FAIL                   = 4, /*!< MPI Comm Duplication failed  */
    COMSIZE_FAIL                   = 5,
    MESH_DIVISIBLE                 = 6,
    MPI_GET_RANK_FAIL              = 7,
    MPI_COMSIZE_FAIL               = 8,
    MPI_DIST_GTAPH_FAIL_X          = 9,
    MPI_DIST_GTAPH_FAIL_Y          = 10,
    GRAPH_CREATE_FAIL              = 11,
    MPI_INEIGHBOR_FAIL_ZX          = 12,
    MPI_INEIGHBOR_FAIL_XY          = 13,
    BLOCK_NUMBER_FAIL              = 14,
    ALLOCATION_FAIL                = 15, /*!< Not enough memory to allocate */
    MPI_FINALIZE_FAIL              = 16,
    MPI_ERROR_HANDLE_FAIL          = 16,
    CONNECTIVITY_CONSTRUCTION_FAIL = 17,
    THOMAS_FAIL                    = 18, /*!< Thomas algorithm failed */
    INPUT_ARGS                     = 19,
    CUFFT_FAIL_X                   = 20, /*!< CUFFT Failure in X-direction transoformation */
    CUFFT_FAIL_INV_X               = 21, /*!< CUFFT Failure in inverse transform in X-direction */
    CUFFT_FAIL_Y                   = 22, /*!< CUFFT Failure in Y-direction transoformation */
    CUFFT_FAIL_INV_Y               = 23, /*!< CUFFT Failure in inverse transformation in Y-direction */

} PittPackResult;

const double pi = 3.1415926535897932384;

// These are for debugging and turning on and off several stages of the Poisson solve

#define POSS 1
#define SOLVE 1
#define IFFTX 1
#define IFFTY 1 
#define FFTX 1
#define FFTY 1
#define USE_SHARED 0 /*!< set to 0 will use global memory, set to 1 will use shared memory*/
#define REVTRSP 1    /*!< 1) stands for simplest transform. 0) transposes using shared mem*/
#define R_COPY 1
#define THOM_FULL_BATCH 0 /*! full batch=1 refers to 2D batch of the system of tridiagonal solves, 0 implies x-dimenional batching only*/
#define JIC 0

PittPackResult OPENACC_Init( int &my_rank, int &com_size );       /*!< Initializes the GPU's. equivalent to MPI_Init()  */
PittPackResult HostToDeviceAssign( int &my_rank, int &com_size ); /*!< Binds each CPU to a GPU  */
const char *   PittPackGetErrorEnum( PittPackResult error );      /*!<Information in text form for Exit Codes */

#endif
