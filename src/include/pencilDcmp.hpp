#ifndef _PENCILDCMP_H_
#define _PENCILDCMP_H_
#include "chunkedArray.hpp"
#include "definitions.h"
#include "mpi.h"
#include "phdf5.hpp"
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <stdlib.h>
#if ( PITTPACKACC )
#include "cufft.h"
#endif
#include "multiGrid.hpp"
#include "signalProc.hpp"
#include "triDiag.hpp"
#include <memory>


/*!    \class PencilDcmp
 *     \brief  This Class performs pencil decomposition and encapsulates all the necessary communications required for FFT transformation
 *
 *
 *
      \f{eqnarray*}{
        \Delta u &=& f  \f}
 *   NOTE: inlining the functions will cause the inherited class not to see the function and hence be undefined in CPU version
 * */

typedef unsigned int uint;

typedef struct shuffle /*!< \brief Struct used to generate iax, iay, jax and jay arrays  */
{
    sint         id;
    vector<sint> nbr;
} shuffle0;

class PencilDcmp
{
    protected:
    SignalProc Sig;
    TriDiag    T;
    MultiGrid  MG;                        /*!< Multigrid object to solve the real part of the transformation */
    MultiGrid  MGC;                       /*! Multigrid object to solve the imaginary part of the transformation */
    int        nx, ny, nz;                /*!<Number of points for each processor (includes all the chunks) */
    int        nxChunk, nyChunk, nzChunk; /*!< dimesnion of each chunk  */
    // since depending on the transformation we will be permutating
    int      p0, p1;          /*<!Number of processors in x and y-directions  */
    int      myRank, comSize; /*!< stores the rank and communicator size */
    int      chunkSize;       /*!< stores chunksize=nxChunk*nyChunk*nzChunk  */
    int      nChunk;          /*!< number of chunks for each processor  */
    int      fullSize;
    int      returnVal;
    int      gangTri;
    double mean;
    MPI_Comm Comm;
    MPI_Comm nbrComm[2]; /*!< stores MPI communicators to be used for neighborhood collectives in x- and y-directions */
    MPI_Comm nbrComm0;   /*!< MPI communicator to use for neighborhood collectives in x-dir */
    MPI_Comm nbrComm1;   /*!< MPI communicator to use for neighborhood collectives in y-dir */
    int *__restrict__ indices        = NULL;
    double *__restrict__ Xbox        = NULL;
    double *__restrict__ coords      = NULL;
    double *tmp                      = NULL;
    double *__restrict__ tmpX        = NULL;
    double *__restrict__ tmpY        = NULL;
    double *__restrict__ tmpMGReal   = NULL;
    double *__restrict__ tmpMGImag   = NULL;
    double *__restrict__ x1          = NULL;
    double *__restrict__ x2          = NULL;
    double *__restrict__ x3          = NULL;
    double *__restrict__ crpcr_lower = NULL;
    double *__restrict__ crpcr_upper = NULL;
    double *__restrict__ crpcr_rhs   = NULL;
    MPI_Request *send_request        = NULL;
    MPI_Status * send_status         = NULL;
    MPI_Request *recv_request        = NULL;
    MPI_Status * recv_status         = NULL;
    double scale; /*!< store the scaling factor for the two inverse transformations depending on the boundary conditons */
    sint * faceTag = NULL;
    double finalErr;
    double *__restrict__ num   = NULL;
    double *__restrict__ denum = NULL;
    double *omega              = NULL; /*!< frequency of the Manufactured Solution*/
    double *__restrict__ subDiag
    = NULL; /*!< allocated with the size of 3, stores only the sub-diagonal elements of the tridiagonal system  */
    double *__restrict__ supDiag
    = NULL; /*!< allocated with the size of 3, stores only the super-diagonal elements of the tridiagonal system  */
    int * length = NULL;
    char *bc     = NULL;              /*!< boundary conditions of the entire domain specified at Xmin,Xmax, Ymin,Ymax, Zmin and Zmax  */
    char  transform[3];               /*!< "S" stands for DST, "C" stands for DCT, "F" stands for Fourier Transform, both "S" or "D", "C"*/
    int *__restrict__ tags    = NULL; /*!< 0 : sine tranform, 1: coisne, 2: Periodic and -1 undefined */
    double *__restrict__ dxyz = NULL; /*!< dx=dxyz[0]; dy=dxyz[1]; dz=dxyz[2]   */
    double *     freqs;
    ChunkedArray P;    /*!< an instance of the chunked array which is the main pointer to hold and perform solve  */
    ChunkedArray R;    /*!< Recieve buff */
    ChunkedArray temp; /*!< temporary array for transposing */
                       //   int *nbrs;
                       //   int *nbrsXY;
    double runTime;
    char   nameAppendix[80];
    sint iaxSize;                  /*!< size of the iax array, we use Compressed Row Storage format to compactly store these indices */
    sint *__restrict__ iax = NULL; /*!< saves the extent for the array*/
    sint *__restrict__ jax = NULL; /*!< saves tha shuffling sequence from local to global for a given index*/
    sint jaxSize;                  /*!< saves tha size of the jax array for copy operations*/
    sint iaySize;                  /*!<size of the iax array, we use Compressed Row Storage format to compactly store these indices */
    sint *__restrict__ iay = NULL; /*!< saves the extent for the array in y-direction*/
    sint *__restrict__ jay = NULL; /*!< saves the shuffling sequence from local to global in y direction*/
    sint jaySize;                  /*!< saves tha size of the jax array for copy operations*/

    int *Nbrs[2]            = {NULL, NULL}; /*!< Nbrs[0] is the stride 1 neighbor where as Nbrs[1] is the stride p0 */
    double *__restrict__ dl = NULL;         /*!< subDiagonal Values of the coefficient matrix */
    double *__restrict__ ds = NULL;         /*!< Diagonal Matrices of the coefficient matrix  */
    double *__restrict__ du = NULL;         /*!<  superDiagonal Matrices of the coefficient matrix  */

  
  void debug();


private:
  void           runInfo();
  void changeLocation();
#if ( PITTPACKACC )
#pragma acc routine vector
#endif
    void saveToDest( const sint source, const sint dest, sint dir );

#if ( PITTPACKACC )
#pragma acc routine vector
#endif
    void saveTmpToDest( const double *__restrict__ tmp, const sint dest, sint dir );

#if ( PITTPACKACC )
#pragma acc routine gang
#endif
    void restoreLocationX(); /*!< restores the data stored at each chunk to its original x-major lay-out  */

#if ( PITTPACKACC )
#pragma acc routine gang
#endif
    void restoreLocationY(); /*!< restores the data stored at each chunk to its original y-major lay-out */

#if ( PITTPACKACC )
#pragma acc routine gang
#endif
    void changeLocationX(); /*!<changes storage layout from local to global for perfoming FFT in x-dir*/

#if ( PITTPACKACC )
#pragma acc routine gang
#endif
    void changeLocationY(); /*!<changes storage layout from local to global for perfoming FFT in y-dir*/

#if ( PITTPACKACC )
#pragma acc routine gang
#endif
    void initializeTrigonometric(); /*!< Initializes the source using \f$sin(\omega x)\f$ and \f$cos(\omega)\f$ bases on the
                                                boundary condition */

    void       changeOwnershipPairwiseExchangeZX(); /*!< Redistributes data by rotating from Z to X direction using pairwise exchange */
    void       changeOwnershipPairwiseExchangeXY(); /*!< Redistributes data by rotating from X to Y direction using pairwise exchange*/
#if ( PITTPACKACC )
#pragma acc routine gang
#endif
    void modifyRhsDirichlet();
#if ( PITTPACKACC )
#pragma acc routine gang
#endif
    void changeLocationXOverlap();

#if ( PITTPACKACC )
#pragma acc routine gang
#endif
    void changeLocationYOverlap();
#if ( PITTPACKACC )
#pragma acc routine gang
#endif
    void preprocessSignalAccordinglyReverse( const int direction, const int counter );

#if ( PITTPACKACC )
#pragma acc routine gang
#endif
    void assignTempY2X( const int chunkId );

#if ( PITTPACKACC )
#pragma acc routine gang
#endif
    void assignTempY2Z( const int chunkId );

#if ( PITTPACKACC )
#pragma acc routine gang
#endif
    void assignBackTempZ2Y( const int chunkId );

#if ( PITTPACKACC )
#pragma acc routine gang
#endif
    void assignTempZ2Y( const int chunkId );

#if ( PITTPACKACC )
#pragma acc routine gang
#endif
    void rearrangeX2Y(); /*!< rearranges storage of data  from x-major to y-major format*/

#if ( PITTPACKACC )
#pragma acc routine gang
#endif
    void rearrangeX2YInverse(); /*!< rearranges storage of data  from y-major back to x-major format*/

#if ( PITTPACKACC )
#pragma acc routine gang
#endif
    void rescale(); /*!< Rescales the inverse transformation, dince FFT transforms by FFTW and CuFFT are not normalized */

#if ( PITTPACKACC )
#pragma acc routine vector
#endif
    int thomas( int i, int j, int dir, int index ); /*!< performs direct solve using Thomas algorithm*/

#if ( PITTPACKACC )
#pragma acc routine gang
#endif
    int solveThm( const int index ); /*!< Performs the solve operaton by calling thomas method */

#if ( PITTPACKACC )
#pragma acc routine gang
#endif
    void solveMG(); /*!< index=0 solves for the real part and index=1 solves for the immaginary part */

#if ( PITTPACKACC )
#pragma acc routine gang
#endif
    void solveMGC(); /*!< index=0 solves for the real part and index=1 solves for the immaginary part */

    void eigenVal( ofstream &myfile ); /*!< prints the eigenvcalues to a file for debugging*/

#if ( PITTPACKACC )
#pragma acc routine
#endif
    PittPackReal getEigenVal( int i, int j ); /*!< Calculates the eigen vale based on the offset, note that, at the tridiagonal solve, j-
                                        direction is the dominant direction if nxChunk==nyChunk && dx=dy, the eigenvalues are symmetric  */
#if ( PITTPACKACC )
#pragma acc routine
#endif
    void thomasPeriodic( int i, int j, int dir,
                         int index ); /*!< Sherman-Morrison version of the Thomas algorithm (also called cyclic Thomas)*/
#if ( PITTPACKACC )
#pragma acc routine
#endif
    void thomasSingleBlock( int i, int j, int dir, int index ); /*!< Singke Block solver used for debugging */
 
#if ( PITTPACKACC )
#pragma acc routine gang
#endif
    void rearrangeY2Z();

#if ( PITTPACKACC )
#pragma acc routine gang
#endif
    void rearrangeY2ZInverse();

#if ( PITTPACKACC )
#pragma acc routine gang
#endif
    void assignBackTempY2Z( const int chunkId );

#if ( PITTPACKACC )
#pragma acc routine gang
#endif
    void postprocessSignalAccordinglyReverse( const int direction, const int counter );

#if ( PITTPACKACC )
#pragma acc routine gang
#endif
    void preprocessSignalAccordingly( const int direction, const int counter );
#if ( PITTPACKACC )
#pragma acc routine gang
#endif
    void postprocessSignalAccordingly( const int direction, const int counter );

#if ( PITTPACKACC )
#pragma acc routine gang
#endif
    void fillInArrayContig( const int i, const int j, const int index );

#if ( PITTPACKACC )
#pragma acc routine gang
#endif
    void fillInArrayBack( const int i, const int j, const int index );

    PittPackResult checkInput();

#if ( PITTPACKACC )
#pragma acc routine gang
#endif
    void setDiag( int i, int j );

#if ( PITTPACKACC )
#pragma acc routine gang
#endif
    void copyAnsFromMG( int index );

#if ( PITTPACKACC )
#pragma acc routine gang
#endif
    void copySourceToMG( int index );

#if ( PITTPACKACC )
#pragma acc routine gang
#endif
    int solveThmBatch( const int index );

#if ( PITTPACKACC )
#pragma acc routine worker
#endif
    void fillInArrayContig( const int i, const int j, int index, double *container );

#if ( PITTPACKACC )
#pragma acc routine worker
#endif
    void fillInArrayBack( const int i, const int j, const int index, double *container );

#if ( PITTPACKACC )
#pragma acc routine seq
#endif
    void thomasLowMem( double *tmpMG, double *rh, double diag, int index );

#if ( PITTPACKACC )
#pragma acc routine vector
#endif
    void pcr( int n, double *a, double *c, double *d ); /*!< Parallel Cyclic Reduction */

#if ( PITTPACKACC )
#pragma acc routine worker
#endif
    void fillInArrayContigNormalize( const int i, const int j, int index, double *container, double eig );

#if ( PITTPACKACC )
#pragma acc routine gang
#endif
    void solveCRP( const int index );

#if ( PITTPACKACC )
#pragma acc routine gang
#endif
    void solvePCR( const int index );

#if ( PITTPACKACC )
#pragma acc routine vector
#endif
    void offDiagNormalize( const int i, const int j, double *lower, double *upper, double eig );

#if ( PITTPACKACC )
#pragma acc routine worker
#endif
    void fillInArrayBack( const int i, const int j, double *container, const int index );

    int createNodalCommunicator();

#if ( PITTPACKACC )
#pragma acc routine seq
#endif
    void imposeBoundaryonContainer( int i, int j, int index, double eig, double *container );

#if ( PITTPACKACC )
#pragma acc routine seq
#endif
    void imposeBoundaryonOffDiag( double eig, double *lower, double *upper );

#if ( PITTPACKACC )
#pragma acc routine seq
#endif
    void imposeBoundaryonCRPTmp( int i, int j, int index, double eig, double *container );

#if ( PITTPACKACC )
#pragma acc routine gang
#endif
    void clear( double *x );
    void print();
    void subtractMeanValue( );
  

  private:
    int        getPeriodicRank( int rank );
    int        getPeriodicIndex( int rank );
    void           getChunkSize();
    void           allocateChunks(); /*!< Allocates the chunkedArray container to hold data */
    void           getNoOfChunk();
    void           MPIStartUp();
    static int getPeriodicRankYZ( int rank );
    void setDxyz( double *X );
    int  getIndex( int rank );
    int  getPeriodicIndexStride( int index );
    int  getPeriodicRankStride( int index );
    void adjustIndex( int &index ); /*!< Adjusts index for the cases where index<0 or  index>p0 */
    int  getTransposeRank();
    void reshape();
    sint oldIndexX( const sint chunkId, const sint j, const sint k );
    sint oldIndexY( const sint chunkId, const sint i, const sint k );
    void decomposeOldIndexX( const sint id, sint &chunkId, sint &j, sint &k );
    void decomposeOldIndexY( const sint id, sint &chunkId, sint &i, sint &k );
    //  sint getDestinationdIndex(const sint id);
    sint getDestinationLocX( const sint id );
    sint getDestinationLocY( const sint id );
    void setUpShuffleArraysX( vector<shuffle0> &a ); /*!< This is only done once on CPU, it sets up the order in which data shuffling should
                                                        be performed in x-dir */

    void nbrAllToAllZX();         /*!< Redistributes data by rotating from Z to X by using neighborhood collective */
    void nbrAllToAllXY();         /*!< Redistributes data by rotating from X to Y by using neighborhood collective */
    void checkGraph( int index ); /*!<Checks the consistency of the graph */

    void constructShuffleVectorX(); /*!<  constructs iay, jay,(Compressed Row Storage) to be copied to GPU for x-major direction*/
    void constructShuffleVectorY(); /*!< constructs iay, jay,(Compressed Row Storage) to be copied to GPU for y-major direction*/
    void setUpShuffleArraysY( vector<shuffle0> &a );
    void saveToTmp( const sint id, double *__restrict__ tmp, sint dir ); /*!< saves the data to a temporary array */
     /*!< This is only done once on CPU, it sets up the order in which data shuffling should be performed in y-dir */
    void restoreTmp( const sint id, double *__restrict__ tmp,
                     sint       dir );     /*!< Assignes the data in the temporary array to its new location */
  
    void printX( ofstream &myfile ); /*!<Each rank prints x-major their own file name data* where '*' is appended by processor rank  */
    void extractTag();                     /*!< Extracts transformation tags based on the Boundary Conditions */

    void printY( ofstream &myfile ); /*!<Each rank prints y-major their own file name data* where '*' is appended by processor rank  */
    void setEigenVal(); /*!< sets the eigenvalue based for each filament  */
    void assignTempX2Yv1( const int chunkId );
    void assignBackTempX2Y( const int chunkId );
    void assignBackTempY2X( const int chunkId );
    void setScale(); /*!< sets the scale of the transformation according to the BC, since FFTW and CUDFFT both do not normalize the
                        tranformation, division by the \f$scale$\f is needed to retrieve the original signal */

    void currentDateTime();
    void detectTransforms(); /*!< Detects the appropriate transformation, i.e. DST00, DCT10, etc based on the boundary conditions */
    void assignTempX2Y( const int chunkId, const int k, double *tmp );

    void nbrAllToAllZXOverlap();
    void nbrAllToAllXYOverlap();

    virtual void performTransformXdir()=0;                /*!< Performs FFT complex transformation in X direction*/
    virtual void performInverseTransformXdir()=0;         /*!< Performs IFFT complex transformation in X direction*/
    virtual void performTransformYdir()=0;                /*!<Perform FFT complex transformation in Y-direction */
    virtual void performInverseTransformYdir()=0;         /*!<Perform IFFT complex transformation in Y-direction */
  
// public interface
    public:
    PencilDcmp(){};                                                 /*!< Class constructor */
    PencilDcmp( int nx, int ny, int nz, int p0, int p1 );           /*!< Alternative Class constructor */
    PencilDcmp( int argcs, char *pArgs[], int nx, int ny, int nz ); /*!< Alternative Class constructor */
    void           setCoords( int dir );    /*!< sets teh coordinates of every pencil depending on the direction of the pencil */
    PittPackResult constructConnectivity(); /*!< neighborhood connectivity to be used in the MPI_Neighbor collectives */
    void graphCreate();           /*!< Creates Neighborhood Graph for communicaton */
    void assignBoundary( char *boundary ); /*!< assig the boundary from user to the object*/
    void assignRhs( double *rhs );
    void fillTrigonometric( double *rhs );
    void readIn();
    void pittPack(); /*!<  method to call FFT and Thomas for solve  */

 
#if ( PITTPACKACC )
#pragma acc routine gang
#endif
    double getError();
    void setBox( double *X );
    void IO( int app, int dir, int aligndir ); /*!<I/O is performed using hdf5 and xmf file formats */
    friend class VTK;                          /*!< for small problems VTK is also provided */

        // rearrnages in z- direction
    virtual ~PencilDcmp(); /*!< Class destructor*/
};

/*!    \class PoissonCPU
 * \brief  This Class inherits Classs PencilDcmp methods and implements the fourier transforms using fftw3
 *
 */

class PoissonCPU : public PencilDcmp
{
public:
    PoissonCPU( int argcs, char *pArgs[], int nx, int ny, int nz ): PencilDcmp( argcs, pArgs, nx, ny, nz ){}; /*!< Alternative Class constructor */
    PoissonCPU( int nx, int ny, int nz, int p0 ) : PencilDcmp( nx, ny, nz, p0, p0 ){}; /*!< constructor */

private:
   void readYLine( int j, fftw_complex *out );
   void writeXLine( int j, fftw_complex *outC ); /*!< reads the appropriate line ans assigns it to fftw for FFT in y-direction*/
   void writeYLine( int j, fftw_complex *outC );
   void readXLine( int j, fftw_complex *out ); /*!< reads the appropriate line ans assigns it to fftw for FFT in x-direction */

    void performTransformXdir() override ;                /*!< Performs FFT complex transformation in X direction*/
    void performInverseTransformXdir() override;         /*!< Performs IFFT complex transformation in X direction*/
    void performTransformYdir() override;                /*!<Perform FFT complex transformation in Y-direction */
    void performInverseTransformYdir() override;         /*!<Perform IFFT complex transformation in Y-direction */
//    void testDST10();
//    void testDST01();

};

/*!    \class PoissonGPU
 * \brief  This Class inherits methods from Classs PencilDcmp and implements the Fourier transforms using NVIDIA's cufft library,
 *
 */

class PoissonGPU : public PencilDcmp
{

public:
    PoissonGPU( int nx, int ny, int nz, int p0 ) : PencilDcmp( nx, ny, nz, p0, p0 ){}; /*!< constructor */
    PoissonGPU( int argcs, char *pArgs[], int nx, int ny, int nz )
    : PencilDcmp( argcs, pArgs, nx, ny, nz ){}; /*!< Alternative Class constructor */
    //    PoissonGPU( int nx, int ny, int nz, int p0 );  /*!< constructor */
    PittPackResult initializeAndBind(); /*was giving strange errors, tried to encapsulate with the constructor but it did not*/

private:
//  the only place we need to override is when we use cufftw
    void performTransformXdir() override;        /*!< Performs FFT complex transformation in X direction*/
    void performInverseTransformXdir() override; /*!< Performs IFFT complex transformation in X direction*/
    void performTransformYdir() override;        /*!<Perform FFT complex transformation in Y-direction */
    void performInverseTransformYdir() override; /*!<Perform IFFT complex transformation in Y-direction */
    void triDiagCusparse( double *dl, double *ds, double *du, double *rhs );
    void debug1();
    // void thomasCusparse(double *rhs);
    // void triDiagCusparse(double *rhs);

};

std::unique_ptr<PencilDcmp>make_Poisson(int argcs, char *pArgs[], int nx, int ny, int nz );/*! factory to make poisson class */

#endif
