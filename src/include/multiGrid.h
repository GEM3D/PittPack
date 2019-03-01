/************************************************************************
 *
 *                           Developed by:
 *
 *                 Jaber Hasbestan and Inanc Senocak
 *
 *
 ************************************************************************
 */

/*!    \class  MultiGrid
 *     \brief  OpenACC implementation of Multigrid algorithm for the solution of the Poisson equation in one dimension
 *
 *
      \f{eqnarray*}{
        \Delta u &=& f  \f}
 *
 *
 * */

#ifndef _MULTIGRID_
#define _MULTIGRID_
#include "math.h"
#include "params.h"
#include "stdio.h"
//#include <stdint.h>
#include <stdlib.h>
#include "chunkedArray.h"
#include "definitions.h"


class MultiGrid
{
    private:
    int     N;
    double *delx;
    int     maxLevel;                    /*!< Maximum level */
    int     arraySize;                   /*!< array size */
    int     tmpArraySize;                /*!< This is only used if Weigted-Jacobi is used as smoother */
    double *__restrict__ u       = NULL; /*!< array containing unkonwns  */
    double *__restrict__ rhs     = NULL; /*!< array containing rhs (or source terms)  */
    double *__restrict__ rhsImag = NULL; /*!< array containing rhs (or source terms)  */
    double *__restrict__ res     = NULL; /*!< array containing residual  */
    double *__restrict__ utmp    = NULL; /*!< array containing residual  */
    double *__restrict__ subDiag0 = NULL; /*!< subDiagonal Values of the coefficient matrix */
    double *__restrict__ supDiag0 = NULL; /*!< superDiagonal Matrices of the coefficient matrix  */
    double *__restrict__ onDiag0  = NULL; /*!<  Diagonal Matrices of the coefficient matrix  */
//    double *__restrict__ dl = NULL; /*!< subDiagonal Values of the coefficient matrix */
//    double *__restrict__ d = NULL; /*!< superDiagonal Matrices of the coefficient matrix  */
//    double *__restrict__ du  = NULL; /*!<  Diagonal Matrices of the coefficient matrix  */


    int *indices
    = NULL; /*!< Save the indices to an array instead of calculating them  each time since there is a recursive sequential loop in the
function*/

    int selectSmoother; /*!<selects smoother 0 for Jacobi and 1 for RB-GS */
    int iterDescend;    /*!< Number of Iterations at smoothing for the coarsening phase of the full V-cycle */
    int iterClimb;      /*!< Number of iterations in smoothing for moving to the denser grid  */
    int iterOuter;      /*!< Number of outer iterations */
    int gridSizePowerTwoPlusOne;        /*!< If the grid size is not divisible, then an MultiGrid is one, it reverts back to Mono Grid */

    public:
    MultiGrid( int n, double delx, char *sm, int iter1, int iter2, int iter3 ); /*!< Class constructor */
    MultiGrid( int n, double delx, int iter1, int iter2, int iter3 );           /*!< Class constructor uses defaul smoother RBGS */
    MultiGrid(){};                                                              /*!< Class constructor uses defaul smoother RBGS */
    void construct( int n, int iter1, int iter2, int iter3, double *sub, double *sup );
    void setDelx( double dx );
    void fillIndices();             /*!< fills out the required indices given a level*/
    void initialize( double *uIn ); /*!< initilaized the internal arrays */
    int  getIndex( int level );     /*!< Calculates the index given a level */
    int  detectMaxLevel( int n );   /*!< calculates the maximum level assuming that the grid size is equal to 2^n+1*/
    void getRequiredSize();         /*!< calculates the required arraySize for the full V-cycle multi-grid operations  */
#pragma acc routine gang
    void multiGridProlong( int level ); /*! prolongs the data back to the denser grid */
#pragma acc routine gang
    void multiGridRestrict( int level ); /*!< restructs residual to the coarser grid */
#pragma acc routine gang
    void residual( int level ); /*!<calculates the residual*/
#pragma acc routine gang
    void smooth( int level, int iter ); /*!< smooth with either Jacobi of Red-Black Gauss-Seidel  */
#pragma acc routine gang
    void solveMono(double onDiag ); /*!< solves using one grid */
#pragma acc routine vector
    double l1Norm(); /*!< Calculates the l1-norm */
#pragma acc routine vector
    double l2Norm(); /*!< Calculates the l2-norm*/
    void   print( double *u, double delx, int level );
    void   print( double delx, int level );
    void   setSmoother( char *sm );
#pragma acc routine gang
    void redBlackMono(double onDiag); /*!< solves using a weighted Jacobi on one grid, this is just for comparison*/
#pragma acc routine gang
    void weightedJacobiMono(); /*!< solves using a weighted Jacobi on one grid, this is just for comparison*/
#pragma acc routine gang
    void weightedJacobiSmoother( int index0, int upperLimit, double d ); /*!< Weighted jacobi is used as smoothe */
#pragma acc routine gang
    void solveMulti( double *out ); /*!< solves using multi grid  */
#pragma acc routine seq
    void exact(); /*!< solve exactly using GS, be careful not to use weighted jacobi with one iteration at the coaresest level */
#pragma acc routine gang
    void redBlackSmoother( int index0, int upperLimit, double d ); /*!< Red Black Gauss Seidel smoother */

    void setIter( int iter1, int iter2, int iter3 ); /*!< sets the internal variables for iteration in the class*/

#pragma acc routine gang
    void reInit();

#pragma acc routine gang
    void assign( double *uIn );

#pragma acc routine 
    void setDiag( double diag );

#pragma acc routine gang
    void redBlackSmoother( int index0, int upperLimit, double d, int level ); /*!< Red Black Gauss Seidel smoother */
#pragma acc routine gang
    void reset();


#pragma acc routine gang
    void reFill( double *out );

//    void   thomas(); /*!< solves using the non-parallelizable Thomas algorithm, will not scale this is just for demonstrating the
//    difference
//                        between the multiGrid and Thomas*/
#pragma acc routine seq
    void thomasLowMem( double *out, double diag,
                       int index ); /*!< For poisson equation uses only the coefiicients, does not save repetitive values in a array*/

#pragma acc routine vector
    void thomasPutBack( double *tmpMG, int index );

//    void thomasCusparse(double *rhs);

#pragma acc routine gang
void fillInArrayContig( ChunkedArray &P,int nChunk,  int nzChunk  , int i,  int j, int index );

#pragma acc routine gang
void fillInArrayBack(ChunkedArray &P, int nChunk,  int nzChunk ,  int i,  int j,  int index );

    friend class PencilDcmp;
    void suggestSize(int a);


    ~MultiGrid();
};

#endif
