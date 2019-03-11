#ifndef _CHUNKEDARRAY_H_
#define _CHUNKEDARRAY_H_
#include "definitions.h"
#include <assert.h>
#include <iostream>
/*!
 *
 * \class ChunkedArray
 * \brief Stores complex number in a a contiguous way and overloads some operators to provide user-friendly interface
 *
 *
 *
 * */

class ChunkedArray
{
    private:
    int    nChunk;                 /*!< number of chunks owned by each processor*/
    int    Nx;                     /*!< Nx is the number of grid point in X for the whole grid for each process  */
    int    Ny;                     /*!< Nx is the number of grid point in Y for the whole grid for each process */
    int    Nz;                     /*!< Nx is the number of grid point in Z for the whole grid for each process */
    int    nx;                     /*!< nx for each chunk, depending on the direction of chunking this might be nx or nx/nChunk  */
    int    ny;                     /*!< ny for each chunk, depending on the direction of chunking this might be ny or ny/nChunk  */
    int    nz;                     /*!< nz for each chunk, depending on the direction of chunking this might be nz or nz/nChunk  */
    int    arraySize;              /*!<Total size of array */
    double Xa, Xb, Ya, Yb, Za, Zb; /*!< Block coordinates to be passed on to Class Phdf5 for visualization*/
    double *__restrict__ P = NULL; /*!< Restricted main pointer to hold data*/
    int    chunkSize       = 0;
    int    orientation;
    double dx, dy, dz;
    complex<double> *__restrict__ Pc=NULL;

    public:
    ChunkedArray(){};                  /*!< Constructor  */
    void           moveHostToDevice(); /*!< Copys from GPU to CPU  */
    void           moveDeviceToHost(); /*!< Copys from GPU to CPU  */
    PittPackResult allocate( int *n, int nbl );
    int            size();
    int            getChunkSize();
    void           setDirection( int dir );
    void           setCoords( double *X );
//    void print();
#if ( PITTPACKACC )
#pragma acc routine
#endif
     double &operator()( int i, int j, int k, int index );

#if ( PITTPACKACC )
#pragma acc routine
#endif
    double &operator()( int i, int j, int k );

#if ( PITTPACKACC )
#pragma acc routine
#endif
    double &operator()( int chunkId, int dir, int i, int j, int k, int index );

#if ( PITTPACKACC )
#pragma acc routine
#endif
    double &operator()( int i, int j, int k, int dir, int index );

#if ( PITTPACKACC )
#pragma acc routine
#endif
    double &operator()( int i );

    //    void moveDeviceToHost(int chunkId);
    //    void moveHostToDevice(int chunkId); /*!< Copys from GPU to CPU  */

    friend class PencilDcmp;
    friend class PoissonCPU;
    friend class PoissonGPU;
    friend class SignalProc;
    friend class TriDiag;
    friend class MpiCom;

    friend class MultiGrid;

    friend class Phdf5;
    //    friend void poisson();
    //#pragma acc routine
    //    int operator()( int i, int j, int k );
    void getAddress( double *rt ); /*!< used for debugging only */
                                   //    void initializeTrigonometric( int *tag );
                                   //    void setDxyz();
    ~ChunkedArray();               /*!< Destructor of the object*/
};

#endif
