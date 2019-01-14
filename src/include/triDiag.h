#ifndef _TRIDIAG_H_
#define _TRIDIAG_H_
#include "chunkedArray.h"
#include "definitions.h"

/*!
 *
 *
 *    \class TriDiag
 *
 *    \brief This class contains the direct solve methods based on the Thomas's algorithm
 *
 **/

class TriDiag
{
    private:
    double *subDiag = NULL; /*!<  lower diagonal values*/
    double *supDiag = NULL; /*!< upper diagonal values*/
    int     nzChunk;        /*!< number of points at z-direction for every chunk */
    int     nChunk;         /*!< number of chunks*/

    public:
    TriDiag(){};                                                  /*!< class constructor */
    void setElems( int nCh, int nzCh, double *sub, double *sup ); /*!< assigns the private variables after construction */
#pragma acc routine
    void thomas( ChunkedArray &P, double *onDiag, const int i, const int j, const int dir,
                 const int index ); /*!< thomas with Dirichlet and Neuman Boundary conditions */
#pragma acc routine
    void thomasPeriodic( ChunkedArray &P, double *onDiag, const int i, const int j, const int dir,
                         const int index ); /*1<Single Block solve for debugging  */
#pragma acc routine
    void thomasSingleBlock( ChunkedArray &P, double *onDiag, const int i, const int j, const int dir,
                            const int index ); /*!<Sherman Morrison modification of the Thomas algorithm for periodic
                                                  boundary condition in the solve direction */
#pragma acc routine
    int thomasReal( ChunkedArray &P, double *onDiag, const int i, const int j,
                    const int dir ); /*!< thomas with Dirichlet and Neuman Boundary conditions */

#pragma acc routine
    void thomasPeriodicReal( ChunkedArray &P, double *onDiag, int i, int j, int dir );

    ~TriDiag(); /*!< destructor of the class*/
};

#endif
