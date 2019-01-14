#ifndef _PARAMS_H_
#define _PARAMS_H_

#ifndef NXCHUNK1
#define NXCHUNK1 64
#endif
#ifndef NYCHUNK1
#define NYCHUNK1 64
#endif
#ifndef NZCHUNK1
#define NZCHUNK1 64
#endif

#define OPENACC _OPENACC
#define EXACT 0
//#define CHUNKSIZE 1000 /*! for defining chunkSize at openACC in temp[chunksize] */
#define DEBUG2 0 /*! turns on debug for poissonGPU */
#define DEBUG_COMM 0
#define DEBUG1 0       /*!<mainly turn off IO for debugging */
#define DEBUG0 0       /*!<mainly turn off IO for debugging  for solverCPU */
#define DEBUG 0        /*!<mainly turn off IO for debugging  for solverCPU */
#define COMM_PATTERN 0 /*!< two type of comm_patterns are available, (0) Pairwise exchange (1) Neighorhood collective */
#define SHORT                                                                                                                              \
    0 /*!< controls the type for \f$iax, iay, jax, jay\f$ arrays, setting to 1 selects \f$short int\f$  where settingt it to zero will set \
the type as \f$int\f$ */
#define ZSIZE 200 /*! shoud be nzChubk*nChunk, not required if the first 4 letters specifying the boundary are not 'P'  */

#define COMM_ON 1

//#define GPUAWARE 1
//#define GPUAWARE2 1

#define DISABLE_CHECKS 1

// const int CHUNK=5;

// Please note that if your ny*nz*sqrt(COM_SIZE) < 65535
// you can use unsigned short instead of int and save memory by the factor of 2
// considering vlock sizes as multiplis of 16, you should be good to go with short int
// if not go to the file definitions.h and change the typedef there

enum PittPackParams /*!<Parameters to set before compiling  */
{
    I_O                         = 1, /*!< set o zero to disable IO */
    SHIFT                       = 1, /*!< Set to 0 for no shifts and 1 to shift*/
    THOM                        = 1,
    nElem                       = NXCHUNK1, /*!<Number of elements in each direction */
    This                        = 0,        /*!<Junk for debugging, will be removed later */
    GPUAW                       = 1,        /*! chunkwise send/recieve for ZX rotation */
    GPUAW2                      = 1,        /*! chunkwise send/recieve for XY rotation */
    MPI_ERROR_DISABLE           = 0,        /*! if set to 1 it will reset MPI_ERROR_FATAl to MPI_ERROR_RETURN */
    PITT_ABORT                  = 0,
    INCLUDE_ERROE_CAL_IN_TIMING = 0, /*! uses an expensive allreduce function misleading to be included in profiling, takes 8 % of 512M
                                        mesh, suggest truning it off for profiling*/
    MULTIGRID = 1,                   /*!< Uses the multigrid algorithm for solution in the z-direction, if turned-off Thomas is used */
};

typedef enum PittPackErrorCodes {
    SUCCESS                        = 0,
    GPU_INIT_FAILURE               = 1,
    SOLVER_ARRAY_INCONSISTENT      = 2,
    MPI_INIT_FAIL                  = 3,
    MPI_DUP_FAIL                   = 4,
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
    ALLOCATION_FAIL                = 15,
    MPI_FINALIZE_FAIL              = 16,
    MPI_ERROR_HANDLE_FAIL          = 16,
    CONNECTIVITY_CONSTRUCTION_FAIL = 17,
    THOMAS_FAIL                    = 18,
    INPUT_ARGS                     = 19,
} PittPackResult;

const double pi = 3.1415926535897932384;

#define POSS 1
#define SOLVE 1
#define IFFTX 1
#define IFFTY 1
#define FFTX 1
#define FFTY 1
#define REV 1
#define JIC 0

PittPackResult OPENACC_Init( int &my_rank, int &com_size ); /*!< Initializes the GPU's. equivalent to MPI_Init()  */
PittPackResult HostToDeviceAssign( int &my_rank, int &com_size );
const char *   PittPackGetErrorEnum( PittPackResult error ); /*!<Information in tex form for Exit Codes */

#endif
