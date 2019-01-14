#ifndef _MPICOMM_H_
#define _MPICOMM_H_

#include "chunkedArray.h"
#include "mpi.h"
#include "statevector.h"
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <stdlib.h>
#define restrict __restrict__

/*!    \class MpiCom
 *     \brief  This Class performs intra-node communications to carry-out
 *     FFT with the help of GPU
 *
 *
 *
 * */

// vector<bitset<M>> mesh;
typedef unsigned int uint;
class MpiCom
{
    private:
    uint N;
    //    ComplexPressure Pn;
    int nx, ny, nz;
    // since depending onthe transformation we will be permutating
    int      p0, p1;
    int      topology[3];
    int      myRank, comSize;
    int      chunkSize;
    double * P0;
    int      nChunk;
    int      fullSize;
    MPI_Comm Comm;
    MPI_Comm nbrComm[2];
    MPI_Comm nbrComm0;
    MPI_Comm nbrComm1;

    ChunkedArray P;
    ChunkedArray R; /*! Recieve buff */
                    //   int *nbrs;
                    //   int *nbrsXY;

    int *Nbrs[2]; /*! Nbrs[0] is the stride 1 neighbor where as Nbrs[1] is the stride p0 */

    public:
    MpiCom(){};
    MpiCom( int nx, int ny, int nz, int p0, int p1 );
    void getChunkSize();
    void changeOwnership();
    void allocateChunks();
    void getNoOfChunk();
    void initialize();
    void constructConnectivity();
    void changeOwnershipPairwiseExchangeZX();
    int  getPeriodicRank( int rank );
    int  getPeriodicIndex( int rank );

    void changeOwnershipPairwiseExchangeXY();
    void changeOwnershipPairwiseExchangeYZ();

    int getPeriodicRankXY( int rank );
    int getPeriodicRankYZ( int rank );

    // neighborhood collectives

    void graphCreate();
    void nbrAllToAllZX();
    void nbrAllToAllXY();
    void nbrAllToAllYZ();
    void checkGraph( int index );

    ~MpiCom();

    friend class VTK;
};

#endif
