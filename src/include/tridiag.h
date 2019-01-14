#ifndef _HEADER_H_
#define _HEADER_H_
#if ( 0 )
#include <algorithm>
#include <bitset>
#include <cstdint>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <stack>
#include <unordered_map>
#endif
#include <assert.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <vector>
using namespace std;

#define restrict __restrict__

/*!    \class  Tridiag System Solver
 *     \brief  This Class encapsulates the iterative and direct solvers
 *
 *
 *
 *      */

// vector<bitset<M>> mesh;

class Tridiag
{
    private:
    uint N;

    public:
    Tridiag( int n );

    void assign( double sub, double diag, double sup, double *a, double *b, double *c );

    void thomas( double *restrict a, double *restrict b, double *restrict c, double *restrict r,
                 double *restrict u ); /*! Good ole Thomas from Numerical Recipes */

    void shermanMorrisonThomas( double *a, double *b, double *c, const double alpha, const double beta, double *r,
                                double *x ); /*! Cyclic Tridiagonal System from Numerical Recipes */

    void multiply( double *a, double *b, double *c, double *x, double *r ); /*! specialized for tridiagonal systems, only ome loop */

    // this multiplication is specialized to handle symmetric SOR where sub and sup are only multiplied by weight which is w/(w*(2-w))
    //
    void multiplyWeighted( double *a, double *b, double *c, double weight, double *x, double *r );

    void conjugateGradientDiagonalPreconditioner( double *sub, double *dig, double *sup, double *b, double *x0, double tolerance,
                                                  int maxiter );

    void conjugateGradient( double *sub, double *dig, double *sup, double *b, double *x0, double tolerance, int maxiter );

    void symmetricSOR( const double *sub, const double *dig, const double *sup, double omega, double *b, double *x0, double tolerance,
                       int maxiter );

    void conjugateGradientSSORPreconditioner( double *sub, double *dig, double *sup, double weight, double *b, double *x0, double tolerance,
                                              int maxiter );

    void lowMemConjugateGradient( double *sub, double *dig, double *sup, double weight, double *b, double *x0, double tolerance,
                                  int maxiter );

    void lowMemSymmetricSOR( double *sub, double *dig, double *sup, double omega, double *b, double *x0, double tolerance, int maxiter );

    void getDiagSSOR( double *sub, double *dig, double *sup, double weight, double *dig1 );

    void lowMemConjugateGradientDiagonalPreconditioner( double *sub, double *dig, double *sup, double weight, double *b, double *x0,
                                                        double tolerance, int maxiter );

    /*! Note that generally the matrix obtained by SSOR called M, is more diagonally dominant so a simple Jacobi prec supuld work well for
     * this section*/

    void GMRES( double *sub, double *dig, double *sup, double *b, int niter, double *x0, double *result, double tol );

    void cs( double *a1, double *a2, double *c, double *s );

    void norm2( int nnode, int nvar, double *vector, double *norm );

    void backward_sub( int neqn, double **A, double *y, double *x );

    int comp_nums( const int *num1, const int *num2 );

    void normalize( int nnode, int nvar, double *vector, double norm );

    void RGMRES( double *sub, double *dig, double *sup, double *b, int niter, double *x0, double *result, double tol, int restart );

    void debug();

    ~Tridiag(){};
};

#endif
