#ifndef _TRIDIAG_H_
#define _TRIDIAG_H_
#include "chunkedArray.hpp"
#include "definitions.h"
#include "mathFunction.hpp"
#include "params.h"
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
    char *  bc = NULL;

    public:
    TriDiag(){};                                                  /*!< class constructor */
    void setElems( int nCh, int nzCh, double *sub, double *sup ); /*!< assigns the private variables after construction */

#if ( PITTPACKACC )
#pragma acc routine
#endif
    void thomas( ChunkedArray &P, double *onDiag, const int i, const int j, const int dir,
                 const int index ); /*!< thomas with Dirichlet and Neuman Boundary conditions */

#if ( PITTPACKACC )
#pragma acc routine
#endif
    void thomasPeriodic( ChunkedArray &P, double *onDiag, const int i, const int j, const int dir,
                         const int index ); /*1<Single Block solve for debugging  */

#if ( PITTPACKACC )
#pragma acc routine
#endif
    void thomasSingleBlock( ChunkedArray &P, double *onDiag, const int i, const int j, const int dir,
                            const int index ); /*!<Sherman Morrison modification of the Thomas algorithm for periodic
                                                  boundary condition in the solve direction */

#if ( PITTPACKACC )
#pragma acc routine
#endif
    int thomasReal( ChunkedArray &P, double *onDiag, const int i, const int j,
                    const int dir ); /*!< thomas with Dirichlet and Neuman Boundary conditions */

#if ( PITTPACKACC )
#pragma acc routine
#endif
    void thomasPeriodicReal( ChunkedArray &P, double *onDiag, int i, int j, int dir );

#if ( PITTPACKACC )
#pragma acc routine vector
#endif
    void pcr( int n, double *a, double *c, double *d );

#if ( PITTPACKACC )
#pragma acc routine seq
#endif
    void thomasLowMem( double *tmpMG, double *rh, double diag, int index );

#if ( PITTPACKACC )
#pragma acc routine seq
#endif
    void thomasLowMemNoBC( double *tmpMG, double *rh, double *diag, int index );

#if ( PITTPACKACC )
#pragma acc routine seq
#endif
    void thomasLowMemNoBCV1( double *tmpMG, double *rh, double *diag, int index );


#if ( PITTPACKACC )
#pragma acc routine seq
#endif
    void thomasLowMemNoBCV2( double *tmpMG, double *rh, double *diag, int index );


#if ( PITTPACKACC )
#pragma acc routine vector
#endif
    void crp( const int n, double offdiag, double *tmpA, double *tmpC, double *tmpRHS, double *ThmA, double *ThmC, double *ThmB,
              double *gam1, double *rhs );

#if ( PITTPACKACC )
#pragma acc routine seq
#endif
    void thomasLowMem( int N, double *a, double *b, double *c, double *r, double *gam );

#if ( PITTPACKACC )
#pragma acc routine seq
#endif
    void shermanMorrisonThomas( double *tmpMG, double *rh, double *rh1, double diag, const double alpha, const double beta, int index );
#if ( PITTPACKACC )
#pragma acc routine seq
#endif
    void shermanMorrisonThomasV1( double *tmpMG, double *rh, double *rh1, double diag, const double alpha, const double beta, int index );


#if ( PITTPACKACC )
#pragma acc routine seq
#endif
    bool checkRhs( const double *rh, const int N );

    void assignBC( char *BC );
#if ( PITTPACKACC )
#pragma acc routine seq
#endif
    void enforceZeroMean( double *tmpMG, double *rh, double *diag, int index );


    ~TriDiag(); /*!< destructor of the class*/
};

#endif
