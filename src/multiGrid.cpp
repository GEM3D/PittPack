/*
 *
 *    Developed by Jaber Hasbestan, PhD
 *
 *    del u= -f; on a one dimensional grid
 *
 *
 */
#include "multiGrid.hpp"
#include "mathFunction.hpp"
#include "params.h"
#include "stdio.h"
#include <iostream>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#define BC 0
#define CELL 0
#define RESPRINT 0

using namespace std;

MultiGrid::MultiGrid( int n, double Delx, char *sm, int iter1, int iter2, int iter3 )
{
    N        = n;
    maxLevel = detectMaxLevel( n );
    getRequiredSize();

    setSmoother( sm );
    printf( "max=%d size=%d smoother=%d N=%d\n", maxLevel, arraySize, selectSmoother, N );
    // allocate
    // allocate on the device

    res     = new double[N];
    rhs     = new double[arraySize];
    rhsImag = new double[arraySize];
    u       = new double[arraySize];
    indices = new int[maxLevel];
    delx    = new double[1];

    delx[0] = 1.0;

    fillIndices();

    setIter( iter1, iter2, iter3 );

    if ( selectSmoother == 0 )
    {
        tmpArraySize = arraySize;
        utmp         = new double[tmpArraySize];
    }
    else
    {
        utmp         = new double[1];
        tmpArraySize = 1;
    }

#if ( PITTPACKACC )
#pragma acc enter data create( this [0:1] )
#pragma acc update device( this [0:1] )
#pragma acc enter data create( rhs [0:arraySize] )
#pragma acc enter data create( rhsImag [0:arraySize] )
#pragma acc enter data create( u [0:arraySize] )
#pragma acc enter data create( utmp [0:tmpArraySize] )
#pragma acc enter data create( indices [0:maxLevel] )
#pragma acc update device( indices [0:maxLevel] )
#pragma acc enter data create( res [0:N] )
#pragma acc enter data create( delx [0:1] )
#pragma acc update device( delx [0:1] )
#endif
}

void MultiGrid::construct( int n, int iter1, int iter2, int iter3, double *sub, double *sup )
{
    char sm[100] = "RB";
    N            = n;
    maxLevel     = detectMaxLevel( n );
    getRequiredSize();

    setSmoother( sm );
    printf( "max=%d size=%d smoother=%d N=%d\n", maxLevel, arraySize, selectSmoother, N );
    // allocate
    // allocate on the device

    res     = new double[N];
    rhs     = new double[arraySize];
    rhsImag = new double[arraySize];
    u       = new double[arraySize];
    indices = new int[maxLevel];
    delx    = new double[1];

    // the right hand side has already been multiplied by delx
    delx[0] = 1.0;

    subDiag0 = new double[3];
    supDiag0 = new double[3];
    onDiag0  = new double[3];

    for ( int i = 0; i < 3; i++ )
    {
        subDiag0[i] = sub[i];
        supDiag0[i] = sup[i];
    }

    fillIndices();

    setIter( iter1, iter2, iter3 );

    if ( selectSmoother == 0 )
    {
        tmpArraySize = arraySize;
        utmp         = new double[tmpArraySize];
    }
    else
    {
        utmp         = new double[1];
        tmpArraySize = 1;
    }

    suggestSize( n );

#if ( PITTPACKACC )
#pragma acc enter data create( this [0:1] )
#pragma acc update device( this [0:1] )
#pragma acc enter data create( rhs [0:arraySize] )
#pragma acc enter data create( rhsImag [0:arraySize] )
#pragma acc enter data create( u [0:arraySize] )
#pragma acc enter data create( utmp [0:tmpArraySize] )
#pragma acc enter data create( indices [0:maxLevel] )
#pragma acc update device( indices [0:maxLevel] )
#pragma acc enter data create( res [0:N] )
#pragma acc enter data create( delx [0:1] )
#pragma acc enter data create( subDiag0 [0:3] )
#pragma acc enter data create( onDiag0 [0:3] )
#pragma acc enter data create( supDiag0 [0:3] )
#pragma acc update device( subDiag0 [0:3] )
#pragma acc update device( supDiag0 [0:3] )
#pragma acc update device( delx [0:1] )
#endif
}

#if ( PITTPACKACC )
#pragma acc routine
#endif
void MultiGrid::setDiag( double diag )
{
    for ( int i = 0; i < 3; i++ )
    {
        onDiag0[i] = diag;
        //  cout<<" onDiag0  "<<onDiag0[i]<<endl;
        //  cout<<" subDiag0  "<<subDiag0[i]<<endl;
        // cout<<" supDiag0  "<<supDiag0[i]<<endl;
    }
}

MultiGrid::MultiGrid( int n, double Delx, int iter1, int iter2, int iter3 )
{
    N        = n;
    maxLevel = detectMaxLevel( n );
    getRequiredSize();

    char sm[100] = "RB";

    setSmoother( sm );
    printf( "maxxLevel =%d size=%d smoother=%d N=%d\n", maxLevel, arraySize, selectSmoother, N );
    // allocate
    // allocate on the device

    res     = new double[N];
    rhs     = new double[arraySize];
    rhsImag = new double[arraySize];
    u       = new double[arraySize];
    indices = new int[maxLevel];
    delx    = new double[1];
    delx[0] = Delx;
    fillIndices();

    setIter( iter1, iter2, iter3 );

    if ( selectSmoother == 0 )
    {
        tmpArraySize = arraySize;
        utmp         = new double[tmpArraySize];
    }
    else
    {
        utmp         = new double[1];
        tmpArraySize = 1;
    }

#if ( PITTPACKACC )
#pragma acc enter data create( this [0:1] )
#pragma acc update device( this [0:1] )
#pragma acc enter data create( rhs [0:arraySize] )
#pragma acc enter data create( rhsImag [0:arraySize] )
#pragma acc enter data create( u [0:arraySize] )
#pragma acc enter data create( utmp [0:tmpArraySize] )
#pragma acc enter data create( indices [0:maxLevel] )
#pragma acc update device( indices [0:maxLevel] )
#pragma acc enter data create( res [0:N] )
#pragma acc enter data create( delx [0:1] )
#pragma acc update device( delx [0:1] )

#endif
}

void MultiGrid::setIter( int iter1, int iter2, int iter3 )
{
    iterDescend = iter1;
    iterClimb   = iter2;
    iterOuter   = iter3;
}

MultiGrid::~MultiGrid()
{
    // allocate on the device

#if ( PITTPACKACC )
#pragma acc exit data delete ( rhs [0:arraySize] )
#pragma acc exit data delete ( rhsImag [0:arraySize] )
#pragma acc exit data delete ( u [0:arraySize] )
#pragma acc exit data delete ( res [0:N] )
#pragma acc exit data delete ( utmp [0:tmpArraySize] )
#pragma acc exit data delete ( indices [0:maxLevel] )
#pragma acc exit data delete ( subDiag0 [0:3] )
#pragma acc exit data delete ( onDiag0 [0:3] )
#pragma acc exit data delete ( supDiag0 [0:3] )
#pragma acc exit data delete ( delx [0:1] )
#pragma acc exit data delete ( this [0:1] )
#endif

    delete[] res;
    delete[] rhs;
    delete[] rhsImag;
    delete[] u;
    delete[] utmp;
    delete[] indices;
    delete[] onDiag0;
    delete[] subDiag0;
    delete[] supDiag0;
    delete[] delx;
}

void MultiGrid::fillIndices()
{
    for ( int i = 0; i < maxLevel; i++ )
    {
        indices[i] = getIndex( i );
        // prin
        printf( "index[%d] = %d \n", i, indices[i] );
    }
    //#pragma acc update device( indices [0:maxLevel] )
}

void MultiGrid::setSmoother( char *sm )
{
    char c1[100] = "JC";
    char c2[100] = "RB";

    if ( strcmp( sm, c1 ) == 0 )
    {
        selectSmoother = 0;
    }
    else if ( strcmp( sm, c2 ) == 0 )
    {
        selectSmoother = 1;
    }

    printf( "smoother set as %s %d\n", sm, selectSmoother );
}

void MultiGrid::initialize( double *uIn )
{
    memset( rhs, 0, arraySize * sizeof( double ) );
    memset( rhsImag, 0, arraySize * sizeof( double ) );
    memset( u, 0, arraySize * sizeof( double ) );
    memset( res, 0, N * sizeof( double ) );

    for ( int i = 0; i < N; i++ )
    {
        rhs[i] = uIn[i];
    }

#if ( PITTPACKACC )
#pragma acc update device( rhs [0:arraySize] )
#pragma acc update device( u [0:arraySize] )
#pragma acc update device( res [0:N] )
#endif
}

#if ( PITTPACKACC )
#pragma acc routine gang
#endif
void MultiGrid::assign( double *uIn )
{
#if ( PITTPACKACC )
#pragma acc loop gang
#endif
    for ( int i = 0; i < N; i++ )
    {
        rhs[i] = -uIn[i];
    }

    //#pragma acc update device( rhs [0:arraySize] )
}

#if ( CELL == 1 )
void MultiGrid::getRequiredSize()
{
    int size = 0;

    for ( int i = 0; i < maxLevel; i++ )
    {
        size += ( ( 1 << ( i + 1 ) ) + 1 );
    }
    arraySize = size;
    printf( "maxlevel =%d the required size %d  \n ", maxLevel, size );
}
#else

void MultiGrid::getRequiredSize()
{
    int size = 0;

    for ( int i = 0; i < maxLevel; i++ )
    {
        size += ( ( 1 << i + 1 ) - 1 );
    }
    arraySize = size;
    printf( "maxlevel =%d the required size %d  \n ", maxLevel, size );
}

#endif

void MultiGrid::setDelx( double dx )
{
    delx[0] = dx;
#if ( PITTPACKACC )
#pragma acc update device( delx[0] )
#endif
}

#if ( CELL == 1 )
int MultiGrid::detectMaxLevel( int n )
{
    int norg  = n;
    int level = 0;

    gridSizePowerTwoPlusOne = 1;

    while ( n > 1 )
    {
        n = n / 2;
        level++;
    }
    if ( ( ( 1 << level ) + 1 ) != ( norg ) && SOLUTIONMETHOD == 4 )
    {
        // printf( "wrong mesh size for multigrid %d nSize should be 2^n+1 %d\n", ( 1 << level ) + 1, norg );
        cout << RED << " Wrong mesh size for multigrid"
             << ", nSize should be 2^n+1 = " << ( 1 << level ) + 1 << RESET << endl
             << endl;
        cout << RED << " switching to monoGrid" << RESET << endl;
        gridSizePowerTwoPlusOne = 0;

        //     exit( 1 );
    }

    return ( level );
}
#else

int MultiGrid::detectMaxLevel( int n )
{
    gridSizePowerTwoPlusOne = 1;

    int norg = n;

    n = n + 2;

    int level = 0;

    while ( n > 1 )
    {
        n = n / 2;
        level++;
    }
    if ( ( ( 1 << level ) - 1 ) != ( norg ) && SOLUTIONMETHOD == 4 )
    {
        printf( "wrong mesh size for multigrid %d nSize should be 2^n+1 %d\n", ( 1 << level ) - 1, norg - 2 );
        gridSizePowerTwoPlusOne = 0;
        exit( 1 );
    }

    return ( level );
}

#endif

// note that this reduction works on 2^n since we have an extra element
// at the end we add that value to the final result that is being stored
// on the first element
// our awesome binary tree reduction algorithm to calculate the l2norm of the residual
#if ( 0 )
#if ( PITTPACKACC )
#pragma acc routine gang
#endif
double MultiGrid::l2Norm()
{
    // double val=0.0;
    double val;

#if ( PITTPACKACC )
#pragma acc loop
#endif
    for ( int j = 0; j < N; j++ )
    {
        res[j] = res[j] * res[j];
    }
    int c    = 1;
    int size = ( N - 1 );

// printf("%d \n", maxLevel);
#if ( PITTPACKACC )
#pragma acc loop seq
#endif
    for ( int k = 0; k < maxLevel; k++ )
    {
        // int k=0;
        c = ( 1 << k );
#if ( PITTPACKACC )
#pragma acc loop gang
#endif
        for ( int j = 0; j < ( size >> ( k + 1 ) ); j++ )
        {
            res[2 * j] += res[2 * j + c];
        }
    }

    val = squareRoot( res[0] + res[N - 1] ) / N;

    return ( val );
}

#else

#if ( PITTPACKACC )
#pragma acc routine vector
#endif
double MultiGrid::l2Norm()
{
    double val = 0.0;
#if ( PITTPACKACC )
#pragma acc loop vector
#endif
    for ( int j = 0; j < N; j++ )
    {
        val += ( res[j] * res[j] );
    }
    val = sqrt( val ) / N;
    return ( val );
}

#endif

#if ( CELL == 1 )
#if ( PITTPACKACC )
#pragma acc routine seq
#endif
int MultiGrid::getIndex( int level )
{
    /*
      if (level > maxLevel) {
        printf("levle is too big\n");
        exit(0);
      }
    */
    int index = 0;

    for ( int i = 0; i < level; i++ )
    {
        index = index + ( 1 << ( maxLevel - i ) ) + 1;
    }
    return ( index );
}
#else

int MultiGrid::getIndex( int level )
{
    int index = 0;

    for ( int i = 0; i < level; i++ )
    {
        index = index + ( 1 << ( maxLevel - i ) ) - 1;
    }
    return ( index );
}

#endif

#if ( CELL == 1 )
#if ( PITTPACKACC )
#pragma acc routine gang
#endif
void MultiGrid::multiGridProlong( int level )
{
    // here we start from maxLevel and go down from there

    int levelOut = level - 1;
    int d        = 1 << ( levelOut );
    int d0       = 1 << ( level );

    int upperLimit  = 2 + ( N - 2 ) / d0 - 1;
    int upperLimit0 = 2 + ( N - 2 ) / d - 1;

    int index  = indices[levelOut];
    int index0 = indices[level];

    // printf("\nprolong Maxlevel =%d index=%d index0=%d d=%d d=%d upperLimit=%d\n", n, index, index0 ,d,d0, upperLimit);

    // this needs to be modfied as well since now the zero values are at the face
    /*
        printf(" level =%d\n",level);
        printf("Prolong indexL = %d replaced by index0L = %d  indexR  = %d  indexR = %d\n",index ,index0,
       index+upperLimit0,index0+upperLimit ); printf("upperLimit = %d upperLimit0= %d\n",upperLimit, upperLimit0); fflush(stdout);
    */
    u[index]               = u[index0];
    u[index + upperLimit0] = u[index0 + upperLimit];

#if ( PITTPACKACC )
#pragma acc loop
#endif
    for ( int i = 1; i < upperLimit; i++ )
    {
        // u[i + index] = u[i * d];
        u[i * 2 + index] += u[i + index0];
        //  u[i*2 -1 + index] =0.5*(u[i-1+index0]+ u[i +index0]);
        //    printf("i  =%d  ii= %d, %lf \n", i * d0+index,i+index0, u[i * d0]);
    }

#if ( PITTPACKACC )
#pragma acc loop
#endif
    for ( int i = 1; i <= upperLimit; i++ )
    {
        // u[i + index] = u[i * d];
        u[i * 2 - 1 + index] += 0.5 * ( u[i + index0] + u[i - 1 + index0] );
        //  u[i*2 -1 + index] =0.5*(u[i-1+index0]+ u[i +index0]);
    }
}
#else

#if ( PITTPACKACC )
#pragma acc routine gang
#endif
void MultiGrid::multiGridProlong( int level )
{
    // here we start from maxLevel and go down from there

    int levelOut = level - 1;
    int d        = 1 << ( levelOut );
    int d0       = 1 << ( level );

    int upperLimit = ( N ) / d0 - 1;

    int upperLimit0 = ( N ) / d - 1;

    int index  = indices[levelOut];
    int index0 = indices[level];

    // printf("================== \n");

    u[index] += 0.5 * ( u[index0] );
    u[index + upperLimit0] += 0.5 * ( u[index0 + upperLimit] );

    //    printf("i  =%d  ii= %d  %d\n", i * 2+ index ,i+index0, i+index0-1);
    /*
        printf("left  u[%d]  u[%d] \n", index , index0);
        printf("right u[%d]  u[%d]  \n", index+upperLimit0 , index0+upperLimit);
        printf("uuperlimits  %d  %d  \n", upperLimit , upperLimit0);
        printf("indeices  %d  %d \n", index , index0);
    */
    // injection part

#if ( PITTPACKACC )
#pragma acc loop
#endif
    for ( int i = 0; i <= upperLimit; i++ )
    {
        u[i * 2 + index + 1] += u[i + index0];
        // printf("FIRST LOOP %d u[%d]  u[%d]=  \n",i, i * 2+ index+1,i+index0);
    }

#if ( PITTPACKACC )
#pragma acc loop
#endif
    for ( int i = 1; i <= upperLimit; i++ )
    {
        // u[i + index] = u[i * d];
        u[i * 2 + index] += 0.5 * ( u[i + index0] + u[i - 1 + index0] );
        //  u[i*2 -1 + index] =0.5*(u[i-1+index0]+ u[i +index0]);
        // printf("SECOND LOOP u[%d]  u[%d]  u[%d]\n", i * 2+ index ,i+index0, i+index0-1);
    }

    // printf("================== \n");
}

#endif

#if ( CELL == 1 )
// restircting
#if ( PITTPACKACC )
#pragma acc routine
#endif
void MultiGrid::multiGridRestrict( int level )
{
    int levelOut = level + 1;

    int d = 1 << ( levelOut );

    int upperLimit = 2 + ( N - 2 ) / d - 1;

    int index  = indices[levelOut];
    int index0 = indices[level];

    int d0 = 1 << ( level );

    int upperLimit0 = 2 + ( N - 2 ) / d0 - 1;

    double coeff = 1. / ( 2. * ( d0 + 1. ) );

    //  printf("---------------coeff=%lf d= %d d0= %d \n ",coeff, d, d0);

    //  delx[0]=1.0;

#if ( BC == 0 )

    rhs[index]              = coeff * ( ( 1. + d0 ) * res[0] + d0 * res[d0] ) / ( ( 1 << level ) * delx[0] );
    rhs[index + upperLimit] = coeff * ( ( 1. + d0 ) * res[N - 1] + d0 * res[N - d0 - 1] ) / ( ( 1 << level ) * delx[0] );

#endif
    /*
        rhs[index]              = coeff * ((1+d0)*res[0]+ d0*res[d0])/((1<<level) *delx[0]);
        rhs[index + upperLimit] = coeff * ((1+d0)*res[N-1]+ d0*res[N-d0-1])/((1<<level) *delx[0]);
    */

    //    printf("Restrict indexL = %d replaced by index0L = %d  indexR  = %d  indexR = %d\n",index ,index0,
    //    index+upperLimit,index0+upperLimit0 );

    // index0 already takes into account the difference in level

#if ( PITTPACKACC )
#pragma acc loop
#endif
    for ( int i = 1; i < upperLimit; i++ )
    {
#if ( INJECT )

        rhs[i + index] = ( res[i * d] );
#else

        //  printf(" res restricting %lf %lf %lf d=%d \n", res[i*d0-1], res[i*d0],res[i*d0+1], d0);
        //  printf(" index %d %d %d d=%d \n", i*d0-2, i*d0, i*d0+2, d0);
        rhs[i + index] = ( 0.25 * ( res[i * d + d / 2] + 2. * res[i * d] + res[i * d - d / 2] ) ) / ( ( 1 << level ) * delx[0] );

#endif
        // for debugging I do direct injection to be able to follow up
        // printf("d =%d %lf \n", i * d, u[i * d]);
    }
}
#else

#if ( PITTPACKACC )
#pragma acc routine gang
#endif
void MultiGrid::multiGridRestrict( int level )
{
    int levelOut = level + 1;

    int d = 1 << ( levelOut );

    //    int upperLimit = 2 + ( N - 2 ) / d - 1;
    int upperLimit = ( N ) / d;

    int index = indices[levelOut];

#if ( PITTPACKACC )
#pragma acc loop
#endif
    for ( int i = 1; i < upperLimit + 1; i++ )
    {
//  printf(" res restricting %lf %lf %lf d=%d \n", res[i*d0-1], res[i*d0],res[i*d0+1], d0);
//  printf(" index %d %d %d d=%d \n", i*d0-2, i*d0, i*d0+2, d0);
#if ( INJECT )

        rhs[index + i - 1] = ( res[i * d - 1] ) / ( ( 1 << level ) * delx[0] );
//       rhs[index+i-1] = ( res[i * d - 1]  ) ;
//       printf(" current index = %d  used locs in res = %d %lf \n",index+i-1,i*d-1,res[i*d-1] );
#else

        rhs[index + i - 1]
        = ( 0.25 * ( res[i * d + d / 2 - 1] + 2. * res[i * d - 1] + res[i * d - d / 2 - 1] ) ) / ( ( 1 << level ) * delx[0] );

//  printf(" current index = %d  used locs in res = %d %d %d %lf \n",index+i-1 ,i*d+d/2-1,i*d-1,i*d-d/2-1, (1<<level)*delx );
//  printf(" current index = %d  used locs in res = %lf %lf %lf %lf \n",rhs[index+i-1] ,res[i*d+d/2-1],res[i*d-1],res[i*d-d/2-1],
// (1<<level)*delx );
#endif
        // for debugging I do direct injection to be able to follow up
        // printf("d =%d %lf \n", i * d, u[i * d]);
    }
}

#endif

#if ( CELL == 1 )
#if ( PITTPACKACC )
#pragma acc routine gang
#endif
void MultiGrid::residual( int level )
{
    // int n=N;
    int d = 1 << level;
    //  int n1 = detectMaxLevel(n);
    //  int index0 = getIndex(level);
    int    index0     = indices[level];
    int    upperLimit = 2 + ( N - 2 ) / d - 1;
    double dx0        = delx[0] * (double)d;
    //  double dx0=delx*d*delx*d;
    // printf("**** level=%d %lf delx=  n=%d upperLimit=%d index0=%d \n",level,delx,n,upperLimit,index0);
    //
    // printf("index0???????=%d upperLimit =%d\n",index0, upperLimit);

    double onDiag0L;
    double onDiag0R;

#if ( BC == 0 )
    //   onDiag0L= onDiag0[0]-((1<<(level+1))-1.0);
    //   onDiag0R= onDiag0[2]-((1<<(level+1))-1.0);

    onDiag0L = onDiag0[0] - ( 1.0 );
    onDiag0R = onDiag0[2] - ( 1.0 );

#endif

#if ( BC == 1 )
    onDiag0L = onDiag0[0] + ( 1.0 );
    onDiag0R = onDiag0[2] + ( 1.0 );
#endif

    res[0] = ( rhs[index0] * dx0 + ( supDiag0[1] * u[index0 + 1] + (onDiag0L)*u[index0] ) / dx0 );

    res[upperLimit]
    = ( rhs[index0 + upperLimit] * dx0 + ( subDiag0[1] * u[index0 + upperLimit - 1] + (onDiag0R)*u[index0 + upperLimit] ) / dx0 );

    /*
        res[0] =( rhs[index0] * dx0 *dx0 + ( supDiag0[1]*u[index0 + 1] +(onDiag0L )* u[index0] ) );
        res[upperLimit] = ( rhs[index0+upperLimit] * dx0*dx0 + ( subDiag0[1]*u[index0+ upperLimit - 1] +(onDiag0R ) * u[index0+upperLimit] )
       );
    */

#if ( PITTPACKACC )
#pragma acc loop
#endif
    for ( int i = 1; i < upperLimit; i++ )
    {
        res[i * d] = ( rhs[i + index0] * dx0
                       + ( supDiag0[1] * u[i + index0 + 1] + onDiag0[1] * u[i + index0] + subDiag0[1] * u[index0 + i - 1] ) / dx0 );
    }

    /*
    for(int i=0;i< upperLimit+1;i++)
    {
    printf("  delx=%lf index0=%d  res %lf\n",delx,index0, res[i] );
    }
    */
}

#else

#if ( PITTPACKACC )
#pragma acc routine gang
#endif
void MultiGrid::residual( int level )
{
    // int n=N;
    int d = 1 << level;
    //  int n1 = detectMaxLevel(n);
    //  int index0 = getIndex(level);
    int    index0     = indices[level];
    int    upperLimit = ( N ) / d - 1;
    double dx0        = delx[0] * (double)d;

    int start;
    start   = ( 1 << level ) - 1;
    int end = N - start - 1;

    // printf("**** level=%d %lf delx=  n=%d upperLimit=%d index0=%d d= %d d=%d \n",level,delx,N,upperLimit,index0,(1<<level)-1,end);
    /* before integration
        res[start] = ( rhs[index0] * dx0 + ( u[index0 + 1] - 2. * u[ index0] ) / dx0 );
        res[end] = ( rhs[index0+upperLimit] * dx0 + ( - 2. * u[upperLimit + index0] + u[index0 + upperLimit -1 ] ) / dx0 );
    */

    res[start] = ( rhs[index0] * dx0 + ( u[index0 + 1] + onDiag0[1] * u[index0] ) / dx0 );
    res[end]   = ( rhs[index0 + upperLimit] * dx0 + ( onDiag0[1] * u[upperLimit + index0] + u[index0 + upperLimit - 1] ) / dx0 );

#if ( PITTPACKACC )
#pragma acc loop
#endif
    for ( int i = 1; i < upperLimit; i++ )
    {
        // res[i * d+start] = ( rhs[i + index0] * dx0 + ( u[i + index0 + 1] - 2. * u[i + index0] + u[index0 + i - 1] ) / dx0 );
        res[i * d + start] = ( rhs[i + index0] * dx0 + ( u[i + index0 + 1] + onDiag0[1] * u[i + index0] + u[index0 + i - 1] ) / dx0 );
        //   printf("%d %d %d %d\n",i*d,i+index0,i+index0+1,i+index0-1);
    }
    /*
        for ( int i = 0; i < N; i++ )
        {
          printf("res(%d)= %lf \n",i,res[i]);
        }
    */
}

#endif

#if ( CELL == 1 )
#if ( PITTPACKACC )
#pragma acc routine gang
#endif
void MultiGrid::smooth( int level, int iter )
{
    int n          = N;
    int d          = 1 << level;
    int upperLimit = 2 + ( n - 2 ) / d - 1;

    int index0 = indices[level];

    for ( int j = 0; j < iter; j++ )
    {
        if ( selectSmoother == 0 && level != ( maxLevel - 1 ) )
        {
            weightedJacobiSmoother( index0, upperLimit, d );
        }
        // else if(selectSmoother==1)
        else
        {
#if ( CELL == 1 )
            redBlackSmoother( index0, upperLimit, d, level );
#else
            redBlackSmoother( index0, upperLimit, d );
#endif
        }
    }
}

#else

#if ( PITTPACKACC )
#pragma acc routine gang
#endif
void MultiGrid::smooth( int level, int iter )
{
    int n          = N;
    int d          = 1 << level;
    int upperLimit = ( N ) / d - 1;

    // printf("D====== %d = level= %d \n ",d,level );

    int index0 = indices[level];

    for ( int j = 0; j < iter; j++ )
    {
        // if ( selectSmoother == 0 && level != ( maxLevel - 1 ) )
        //  if ( selectSmoother == 0  )
        {
        //       weightedJacobiSmoother( index0, upperLimit, d );
        } // else if(selectSmoother==1)
        //   else
        {
            redBlackSmoother( index0, upperLimit, d );
        }
    }
}

#endif

#if ( CELL == 1 )
#if ( PITTPACKACC )
#pragma acc routine gang
#endif
void MultiGrid::weightedJacobiSmoother( int index0, int upperLimit, double d )
{
#if ( PITTPACKACC )
#pragma acc loop
#endif
    for ( int i = 0; i < upperLimit + 1; i++ )
    {
        utmp[i] = u[i + index0];
    }

#if ( PITTPACKACC )
#pragma acc loop
#endif
    for ( int i = 1; i < upperLimit; i = i + 1 )
    {
        u[i + index0] = ( rhs[i + index0] * d * delx[0] + ( supDiag0[1] * utmp[i + 1] + subDiag0[1] * utmp[i - 1] ) / ( d * delx[0] ) ) * d
                        * delx[0] * ( -1. / onDiag0[1] ) * 2. / 3.
                        + 1. / 3. * utmp[i];
    }
}

#else

#if ( PITTPACKACC )
#pragma acc routine gang
#endif
void MultiGrid::weightedJacobiSmoother( int index0, int upperLimit, double d )
{
    int farRight = index0 + upperLimit;

// printf(" smoother %d tmpArraySize= %d farRight = %d d=%lf \n",selectSmoother,tmpArraySize,farRight,d);
#if ( PITTPACKACC )
#pragma acc loop
#endif
    /* for ( int i = 0; i < arraySize; i++ )
     {
         utmp[i] = u[i];
         printf("********** utmp[%d] = %lf u[%d]=%lf\n",i,utmp[i],i, u[i]);
     }
 */
    double coeff = -1. / onDiag0[1];

    if ( index0 != arraySize - 1 )
    {
        u[index0] = ( rhs[index0] * delx[0] * delx[0] * d * d + utmp[index0 + 1] ) * (coeff)*2. / 3.0 + 1. / 3.0 * utmp[index0];
        u[index0 + upperLimit] = ( rhs[index0 + upperLimit] * d * d * delx[0] * delx[0] + utmp[index0 + upperLimit - 1] ) * (coeff)*2. / 3.
                                 + 1. / 3. * utmp[index0 + upperLimit];
    }
    // printf("!!!!!!!!!!!!!!!! index0 = %d  upper =%d  d=%lf  u[%d]=%lf \n",index0, index0+upperLimit, d,farRight,utmp[farRight]);

#if ( PITTPACKACC )
#pragma acc loop
#endif
    for ( int i = 1; i < upperLimit; i = i + 1 )
    {
        u[i + index0] = ( rhs[i + index0] * d * delx[0] + ( utmp[index0 + i + 1] + utmp[index0 + i - 1] ) / ( d * delx[0] ) ) * d * delx[0]
                        * (coeff)*2. / 3.
                        + 1. / 3. * utmp[i + index0];
        // printf(" index0 = %d  right =%d  left =%d u=%lf rhs= %lf \n",index0+1,index0+1+i,index0-1+i,
        // upperLimit,u[i+index0],rhs[i+index0]);
    }

    // printf(" index0 = %d  right =%d  left =%d \n",index0,index0+1,index0-1);
}

#endif

#if ( CELL == 1 )

#if ( PITTPACKACC )
#pragma acc routine gang
#endif
void MultiGrid::redBlackSmoother( int index0, int upperLimit, double d, int level )
{
    int    coeff = ( 1 << ( level + 1 ) ) - 1;
    double onDiag0L;
    double onDiag0R;

#if ( BC == 0 )
    /*
       onDiag0L= onDiag0[0]-((1<<(level+1))-1.0);
       onDiag0R= onDiag0[2]-((1<<(level+1))-1.0);

    */
    // for debug
    onDiag0L = onDiag0[0] - ( 1.0 );
    onDiag0R = onDiag0[2] - ( 1.0 );

#elif ( BC == 1 )
    onDiag0L = onDiag0[0] + 1.;
    onDiag0R = onDiag0[2] + 1.;
#endif

    //  printf(" ondiagL =%lf level= %d coeff=%d \n",onDiag0[1],level,coeff );
    //  delx[0]=1.0;

    u[index0] = ( rhs[index0] * d * d * delx[0] * delx[0] + ( supDiag0[1] * u[index0 + 1] ) ) * ( -1. / onDiag0L );
    u[index0 + upperLimit]
    = ( rhs[index0 + upperLimit] * d * d * delx[0] * delx[0] + ( subDiag0[1] * u[index0 + upperLimit - 1] ) ) * ( -1. / onDiag0R );

#if ( PITTPACKACC )
#pragma acc loop
#endif
    for ( int i = 1; i < upperLimit / 2; i++ )
    {
        u[2 * i + index0] = ( rhs[2 * i + index0] * d * d * delx[0] * delx[0]
                              + ( supDiag0[1] * u[2 * i + index0 + 1] + subDiag0[1] * u[2 * i + index0 - 1] ) )
                            * ( -1. / onDiag0[1] );
    }

#if ( PITTPACKACC )
#pragma acc loop
#endif
    for ( int i = 1; i < upperLimit; i = i + 2 )
    {
        u[i + index0]
        = ( rhs[i + index0] * d * d * delx[0] * delx[0] + ( supDiag0[1] * u[i + index0 + 1] + subDiag0[1] * u[i + index0 - 1] ) )
          * ( -1. / onDiag0[1] );
    }
}
#else

#if ( PITTPACKACC )
#pragma acc routine gang
#endif
void MultiGrid::redBlackSmoother( int index0, int upperLimit, double d )
{
    // at each level the end elements are calculated

    double coeff = -1. / onDiag0[1];

#if ( PITTPACKACC )
#pragma acc loop
#endif
    for ( int i = 1; i < upperLimit; i = i + 2 )
    {
        //    uold=u[i+index0];
        u[i + index0] = ( rhs[i + index0] * d * d * delx[0] * delx[0] + ( u[i + index0 + 1] + u[i + index0 - 1] ) ) * coeff;
        // printf("[ i + index0] = %d , [i + index0 +1] = %d ,[i + index0 - 1 ]=%d \n",i + index0, i + index0+1, i + index0-1);
    }
    if ( index0 != arraySize - 1 )
    {
        u[index0] = ( rhs[index0] * delx[0] * delx[0] * d * d + u[index0 + 1] ) * coeff;
    }
#if ( PITTPACKACC )
#pragma acc loop
#endif
    for ( int i = 1; i < upperLimit / 2; i++ )
    {
        u[2 * i + index0] = ( rhs[2 * i + index0] * d * d * delx[0] * delx[0] + ( u[2 * i + index0 + 1] + u[2 * i + index0 - 1] ) ) * coeff;

        // printf("[2 * i + index0] = %d , [2 * i + index0 +1] = %d , [2 * i + index0 - 1 ]=%d \n",2 * i + index0,2 * i + index0+1,2 * i +
        // index0-1);
    }
    if ( index0 != arraySize - 1 )
    {
        u[index0 + upperLimit] = ( rhs[index0 + upperLimit] * d * d * delx[0] * delx[0] + u[index0 + upperLimit - 1] ) * coeff;
    }

    /*
    for(int i=0;i<arraySize;i++)
    {
    printf("u[%d]=%lf\n",i,u[i]);
    }
    */
}

#endif

void MultiGrid::print( double *u, double delx, int level )
{
    int d          = 1 << level;
    int upperLimit = 2 + ( N - 2 ) / d;

    // printf("upper limit =%d \n",upperLimit);
    for ( int i = 0; i < upperLimit; i++ )
    {
        // printf("i= %d %lf \n",i*d,u[i*d]);
        printf( "%lf %lf \n", i * d * delx, u[i * d] );
    }
}

void MultiGrid::print( double delx, int level )
{
    int d          = 1 << level;
    int upperLimit = 2 + ( N - 2 ) / d;

    // printf("upper limit =%d \n",upperLimit);
    for ( int i = 0; i < upperLimit; i++ )
    {
        // printf("i= %d %lf \n",i*d,u[i*d]);
        //    printf("%lf %lf \n", i * d * delx[0], res[i * d]);
        printf( "%lf %lf \n", i * d * delx, u[i * d] );
    }
}

#if ( PITTPACKACC )
#pragma acc routine gang
#endif
void MultiGrid::reFill( double *out )
{
#if ( PITTPACKACC )
#pragma acc loop gang
#endif
    for ( int i = 0; i < N; i++ )
    {
        out[i] = u[i];
        // out[i] = 7.*25.;
        //    cout<< " answer "<< out[i]<<"delx[0] "<<delx[0]<<endl;
    }
}

#if ( PITTPACKACC )
#pragma acc routine gang
#endif
void MultiGrid::solveMono( double onDiag )
{
#if ( 1 )
    int level      = 0;
    int d          = 1 << level;
    int n          = N;
    int upperLimit = N;
    int index0     = 0.0;

    //  printf("start =%d upperLimit= %d \n",index0,upperLimit );
    // printf("\n inside res start =%d upperLimit= %d delx[0]=%lf \n",index0,upperLimit, d*delx[0] );

    double resTot = 1.0;
    int    count  = 0;

    //    printf( "Mono solving" );

    double val = 0.0;
    residual( 0 );

    resTot = l2Norm();

    if ( resTot < 1.e-10 )
    {
        return;
    }

    //            redBlackMono();
    while ( resTot > 1.e-10 )
    //#pragma acc loop seq
    // for(int k=0;k<20;k++)
    {
        // weightedJacobiMono();
        redBlackMono( onDiag );
        residual( level );
        //      count++;
        //     resTot = 0.0;
        resTot = l2Norm();
// val = 0.0;

/*
#pragma acc loop vector reduction( + : val )
        for ( int j = 0; j < N; j++ )
        {
            val += ( res[j] * res[j] );
        }

        resTot = val / N;
*/

//        count++;
#if ( RESPRINT == 1 )

        printf( "monitor res %lf %e\n ", onDiag, resTot );
#endif

        //        if ( count % 100 == 0 )
        //        {
        //    printf("monitor res %e \n",resTot);
        //        }
    }
#endif

#if ( RESPRINT == 1 )
    printf( "count=%d onDiag = %lf \n", count, onDiag );
#endif
    // printf();

    //    reFill( out );
    /*
        #pragma acc loop vector
         for(int i=0;i<N;i++)
        {
             out[i]=10*25;
        //    cout<< " answer "<< out[i]<<"delx[0] "<<delx[0]<<endl;
        }
      */
    /*
        // print(u,delx[0], level);
        if(count>1)
    {
         printf("number of iterations for full blown RB-GS %d %e\n",count, resTot);
    }
    */
}

#if ( CELL == 1 )
#if ( PITTPACKACC )
#pragma acc routine gang
#endif
void MultiGrid::redBlackMono( double onDiag )
{
    int    upperLimit = N;
    double onDiagL;
    double onDiagR;
    int    rem = N % 2;

    if ( rem == 0 )
    {
        rem = -1;
    }

    //   double onDiag0D=onDiag0[1];
    //  cout<< " REDBLACK "<<endl;

#if ( 1 )

#if ( BC == 0 )
    //   onDiagL = onDiag - 1.;
    //   onDiagR = onDiag - 1.;

    u[0] = ( rhs[0] + ( supDiag0[1] * u[1] ) ) * ( -1. / ( onDiag0[1] - 1.0 ) );

    u[N - 1] = ( rhs[N - 1] + ( subDiag0[1] * u[N - 2] ) ) * ( -1. / ( onDiag0[1] - 1.0 ) );

#endif

// Starting from red, this loop handles red cells first
#if ( PITTPACKACC )
#pragma acc loop
#endif
    for ( int i = 1; i < upperLimit / 2; i++ )
    {
        u[2 * i] = ( rhs[2 * i] + ( supDiag0[1] * u[2 * i + 1] + subDiag0[1] * u[2 * i - 1] ) ) * ( -1. / onDiag0[1] );
    }

// revisisted this loop sice PGI had a hard time understanding that we jump one at a time

// Starting from red, this loop handles black cells first
#if ( PITTPACKACC )
#pragma acc loop
#endif
    for ( int i = 1; i < ( upperLimit + rem ); i = i + 2 )
    {
        u[i] = ( rhs[i] + ( supDiag0[1] * u[i + 1] + subDiag0[1] * u[i - 1] ) ) * ( -1. / onDiag0[1] );
    }

#endif

    /*
    #if ( PITTPACKACC )
    #pragma acc loop gang
    #endif
        for ( int i = 0; i < N ; i++ )
        {
            u[i] =-rhs[i];
        }
    */

    /*
     // PGI rejects parallelization of the following although there is no data dependency
     #if ( PITTPACKACC )
     #pragma acc loop
     #endif
         for ( int i = 1; i < upperLimit; i = i + 2 )
         {
             //u[i] = ( rhs[i] * delx[0] * delx[0] + ( supDiag0[1]*u[i + 1] + subDiag0[1]*u[i - 1] ) )*( -1./onDiag0[1]);
             u[i] = ( rhs[i] * delx[0] * delx[0] + ( supDiag0[1]*u[i + 1])+subDiag0[1]*u[i-1])*( -1./onDiag0[1]);
         }
    */
}

#else

#if ( PITTPACKACC )
#pragma acc routine gang
#endif
void MultiGrid::redBlackMono( double onDiag )
{
    int    upperLimit = N;
    double coeff      = -1. / onDiag;

    u[0] = ( rhs[0] * delx[0] * delx[0] + u[1] ) * coeff;

    u[N - 1] = ( rhs[N - 1] * delx[0] * delx[0] + u[N - 2] ) * coeff;

#if ( PITTPACKACC )
#pragma acc loop
#endif
    for ( int i = 1; i < upperLimit / 2; i++ )
    {
        u[2 * i] = ( rhs[2 * i] * delx[0] * delx[0] + ( u[2 * i + 1] + u[2 * i - 1] ) ) * coeff;
    }
#if ( PITTPACKACC )
#pragma acc loop
#endif
    for ( int i = 1; i < upperLimit; i = i + 2 )
    {
        u[i] = ( rhs[i] * delx[0] * delx[0] + ( u[i + 1] + u[i - 1] ) ) * coeff;
    }

    /*
    for(int i=0;i<N;i++)
    {
    printf("u(%d)=%lf\n",i+1,u[i]);
    }
    */
}

#endif

#if ( PITTPACKACC )
#pragma acc routine gang
#endif
void MultiGrid::weightedJacobiMono()
{
#if ( PITTPACKACC )
#pragma acc loop
#endif
    for ( int i = 0; i < N; i++ )
    {
        utmp[i] = u[i];
    }

#if ( PITTPACKACC )
#pragma acc loop
#endif
    for ( int i = 1; i < N - 1; i = i + 1 )
    {
        u[i] = ( rhs[i] * delx[0] * delx[0] + ( utmp[i + 1] + utmp[i - 1] ) ) * 0.5 * 2. / 3. + 1. / 3. * utmp[i];
    }
}

/*
void getErr( double *u, int n, double delx[0], double pi )
{
    double error = 0.0;

    for ( int i = 0; i < n; i++ )
    {
        error += ( u[i] - sin( OMEGA * i * delx[0] * pi ) ) * ( u[i] - sin( OMEGA * i * delx[0] * pi ) );
        // error +=(u[i]-sin(i*delx[0]*pi))*(u[i]-sin(i*delx[0]*pi));
    }


    for(int i=0;i<n;i++)
    {
    if(error<fabs(u[i]-sin(OMEGA*i*delx[0]*pi)))
    {
     error =fabs(u[i]-sin(OMEGA*i*delx[0]*pi));
    }
    }

    printf( "\nError %e delx[0]=%lf \n", error, delx[0] );

}
*/

#if ( CELL == 1 )
#if ( PITTPACKACC )
#pragma acc routine seq
#endif
void MultiGrid::exact()
{
    // smooth( maxLevel - 1, 5 );
    //
    //
    //
    //
    int level = maxLevel - 1;

    int    coeff = ( 1 << ( level + 1 ) ) - 1;
    double onDiag0L;
    double onDiag0R;
    int    n          = N;
    int    d          = 1 << level;
    int    upperLimit = 2 + ( n - 2 ) / d - 1;
    int    index0     = indices[level];
#if ( BC == 0 )
    /*
       onDiag0L= onDiag0[0]-((1<<(level+1))-1.0);
       onDiag0R= onDiag0[2]-((1<<(level+1))-1.0);
    */
    // for debug
    onDiag0L = onDiag0[0] - ( 1.0 );
    onDiag0R = onDiag0[2] - ( 1.0 );

#elif ( BC == 1 )
    onDiag0L = onDiag0[0] + 1.;
    onDiag0R = onDiag0[2] + 1.;
#endif

    //  printf(" ondiagL =%lf level= %d coeff=%d \n",onDiag0[1],level,coeff );
    //  delx[0]=1.0;
    for ( int k = 0; k < 5; k++ )
    {
        u[index0] = ( rhs[index0] * d * d * delx[0] * delx[0] + ( supDiag0[1] * u[index0 + 1] ) ) * ( -1. / onDiag0L );
        u[index0 + upperLimit]
        = ( rhs[index0 + upperLimit] * d * d * delx[0] * delx[0] + ( subDiag0[1] * u[index0 + upperLimit - 1] ) ) * ( -1. / onDiag0R );
        u[index0 + upperLimit - 1]
        = ( rhs[index0 + upperLimit - 1] * d * d * delx[0] * delx[0] + ( supDiag0[1] * u[index0 + upperLimit] + subDiag0[1] * u[index0] ) )
          * ( -1. / onDiag0[1] );
    }
}
#else

#if ( PITTPACKACC )
#pragma acc routine gang
#endif
void MultiGrid::exact()
{
    double dx        = ( ( 1 << ( maxLevel - 1 ) ) * delx[0] );
    u[arraySize - 1] = rhs[arraySize - 1] * dx * dx * ( -1. / onDiag0[1] );
}

#endif

#if ( PITTPACKACC )
#pragma acc routine gang
#endif
void MultiGrid::solveMulti( double *out )
{
    double resTot = 1.0;
    int    count  = 0;

    // printf("Multigrid Solving\n ");
    //   printf("delx inside multigrid= %lf ", delx[0]);

    int bol = 1;

    residual( 0 );

    resTot = l2Norm();

#if ( RESPRINT == 1 )
    printf( "initial res in multigrid= %lf \n", resTot );
#endif

#if ( 1 )

    if ( resTot < 1.e-10 )
    {
        return;
    }

#if ( 1 ) // printf("solving\n");
    while ( resTot > 1.e-10 )
    //         for ( int k = 0; k < 3 ; k++ )
    {
        // decsend
        for ( int i = 0; i < maxLevel - 1; i++ )
        {
            smooth( i, iterDescend );
            residual( i );
            multiGridRestrict( i );
        }

        // this is an issue with weighted jacobi, use regular jacobi for this, solve exactly
        // this section solves exactly

        exact();

        for ( int i = maxLevel - 1; i > 0; i-- )
        {
            smooth( i, iterClimb );
            multiGridProlong( i );
        }

        smooth( 0, iterDescend );

        residual( 0 );

        reInit();

        count++;

        // if(count%2==0)
        {
            resTot = l2Norm();
        }

#if ( RESPRINT == 1 )
        printf( "res=%lf onDiag = %lf \n", resTot, onDiag0[1] );
#endif

        /*
                   if(resTot>1.e10)
                   {
                     printf("resTot=%lf\n",resTot);
                     exit(1);
                    }
        */
        // resTot=l2Norm()/N;
        // resTot=l1Norm();
    }
#endif

#endif
    reFill( out );

#if ( RESPRINT == 1 )
    if ( count > 1 )
    {
        printf( "count in multi=%d onDiag = %lf res=%lf \n", count, onDiag0[1], resTot );
    }
#endif
    /*
    #pragma acc loop vector
    for(int i=0;i<N;i++)
    {
    out[i]=u[i];
    }
    */
    /*
    #if ( PITTPACKACC )
    #pragma acc update self( u [0:arraySize] )
    #endif
    */
    //    print( u, delx[0], 0 );
    //    printf( "Error =%e\n", resTot  );
}

#if ( PITTPACKACC )
#pragma acc routine gang
#endif
void MultiGrid::reInit()
{
#if ( PITTPACKACC )
#pragma acc loop
#endif
    for ( int i = N; i < arraySize; i++ )
    {
        u[i] = 0.0;
    }
}

#if ( PITTPACKACC )
#pragma acc routine gang
#endif
void MultiGrid::reset()
{
#if ( PITTPACKACC )
#pragma acc loop
#endif
    for ( int i = 0; i < arraySize; i++ )
    {
        rhs[i] = 0.0;
        u[i]   = 0.0;
    }

#if ( PITTPACKACC )
#pragma acc loop
#endif
    for ( int i = 0; i < N; i++ )
    {
        res[i] = 0.0;
    }
}
#if ( 0 )
#pragma acc routine seq
void                MultiGrid::thomasLowMem( double *tmpMG )
{
    double bet;
    /*
        double a[3],b[3],c[3];


        for(int i=0;i<3;i++)
        {
          a[i]=subDiag0[i];
          b[i]=onDiag0[i];
          c[i]=supDiag0[i];
        }
    */
    // for Dirirchlet

    //     onDiag0[1]=diag;
    onDiag0[0] = onDiag0[1] - 1.;
    onDiag0[2] = onDiag0[1] - 1.;

    //      b[0]=b[1]-1.;
    //     b[2]=b[1]-1.;
    rhs[0] = rhs[0] / ( bet = onDiag0[0] );
    // cout << r[0] << eNdl;
    // priNtf("r[0]=%lf\N",r[0]);

    int j    = 1;
    tmpMG[j] = supDiag0[j - 1] / bet;
    bet      = onDiag0[1] - subDiag0[1] * tmpMG[j];
    rhs[1]   = ( rhs[1] - subDiag0[1] * rhs[j - 1] ) / bet;

#if ( 1 )
    for ( int j = 2; j < N - 1; j++ )
    {
        //       DecompositioN subDiag0Nd forwsubDiag0rd substitutioN.
        tmpMG[j] = supDiag0[1] / bet;
        bet      = onDiag0[1] - subDiag0[1] * tmpMG[j];
        rhs[j]   = ( rhs[j] - subDiag0[1] * rhs[j - 1] ) / bet;
    }

    j        = N - 1;
    tmpMG[j] = supDiag0[1] / bet;
    bet      = onDiag0[2] - subDiag0[2] * tmpMG[j];
    rhs[j]   = ( rhs[j] - subDiag0[2] * rhs[j - 1] ) / bet;

    //  cout << subDiag0[2] << " " << b[2] << eNdl;

#pragma acc loop seq
    for ( int j = ( N - 2 ); j >= 0; j-- )
    {
        rhs[j] -= tmpMG[j + 1] * rhs[j + 1];
        // cout << " j " << j << eNdl;
    }
#endif
    // printf("thomasing\n");
}
#endif

#if ( PITTPACKACC )
#pragma acc routine seq
#endif
void MultiGrid::thomasLowMem( double *tmpMG, double diag, int index )
{
    double  bet;
    double *rh[2];

    rh[0] = rhs;
    rh[1] = rhsImag;

    double a[3], b[3], c[3];

    int n = N;

    for ( int i = 0; i < 3; i++ )
    {
        a[i] = subDiag0[i];
        //      b[i]=onDiag0[i];
        c[i] = supDiag0[i];
    }

    b[1] = diag;

    // for Dirirchlet

    b[0] = b[1] - 1.;
    b[2] = b[1] - 1.;

    /*
        rh[0] = rh[0] / ( bet = b[0] );
        // cout << r[0] << endl;
        // printf("r[0]=%lf\n",r[0]);

        int j  = 1;
        tmpMG[j] = c[j - 1] / bet;
        bet    = b[1] - a[1] * tmpMG[j];
        rh[1] = ( rh[1] - a[1] * rh[j - 1] ) / bet;

        for ( int j = 2; j < n - 1; j++ )
        {
            //       Decomposition and forward substitution.
            tmpMG[j] = c[1] / bet;
            bet    = b[1] - a[1] * tmpMG[j];
            rh[j] = ( rh[j] - a[1] * rh[j - 1] ) / bet;
        }

        j      = N - 1;
        tmpMG[j] = c[1] / bet;
        bet    = b[2] - a[2] * tmpMG[j];
        rh[j] = ( rh[j] - a[2] * rh[j - 1] ) / bet;

        //  cout << a[2] << " " << b[2] << endl;
        for ( int j = ( n - 2 ); j >= 0; j-- )
        {
            rh[j] -= tmpMG[j + 1] * rh[j + 1];
            // cout << " j " << j << endl;
        }
    */

    rh[index][0] = rh[index][0] / ( bet = b[0] );
    // cout << r[0] << eNdl;
    // priNtf("r[0]=%lf\n",r[0]);

    int j        = 1;
    tmpMG[j]     = c[j - 1] / bet;
    bet          = b[1] - a[1] * tmpMG[j];
    rh[index][1] = ( rh[index][1] - a[1] * rh[index][j - 1] ) / bet;

    for ( int j = 2; j < n - 1; j++ )
    {
        //       DecompositioN and forward substitution.
        tmpMG[j]     = c[1] / bet;
        bet          = b[1] - a[1] * tmpMG[j];
        rh[index][j] = ( rh[index][j] - a[1] * rh[index][j - 1] ) / bet;
    }

    j            = N - 1;
    tmpMG[j]     = c[1] / bet;
    bet          = b[2] - a[2] * tmpMG[j];
    rh[index][j] = ( rh[index][j] - a[2] * rh[index][j - 1] ) / bet;

    //  cout << a[2] << " " << b[2] << eNdl;
    for ( int j = ( n - 2 ); j >= 0; j-- )
    {
        rh[index][j] -= tmpMG[j + 1] * rh[index][j + 1];
        // cout << " j " << j << eNdl;
    }
}

/*
void MultiGrid::thomasCusparse(double *rhs)
{

#pragma acc host_data use_device(rhs, du, dl, d)
  {
    cusparseHandle_t handle = NULL;
    cusparseCreate(&handle);
    cusparseDgtsv_nopivot(handle, N, 1, dl, d, du, rhs, N);
    cusparseDestroy(handle);
  }


}

*/

#if ( PITTPACKACC )
#pragma acc routine vector
#endif
void MultiGrid::thomasPutBack( double *tmpMG, int index )
{
    double *rh[2];

    rh[0] = rhs;
    rh[1] = rhsImag;

#if ( PITTPACKACC )
#pragma acc loop vector
#endif
    for ( int i = 0; i < N; i++ )
    {
        tmpMG[i] = -rh[index][i];
        //   tmpMG[i]=onDiag0[1];
    }
}

#if ( PITTPACKACC )
#pragma acc routine gang
#endif
void MultiGrid::fillInArrayContig( ChunkedArray &P, int nChunk, int nzChunk, int i, int j, int index )
{
//    cout<<"before inside fill in"<<endl;
#if ( PITTPACKACC )
#pragma acc loop gang
#endif
    for ( int id = 0; id < nChunk; id++ )
    {
#if ( PITTPACKACC )
#pragma acc loop vector
#endif
        for ( int k = 0; k < nzChunk; k++ )
        {
            rhs[id * nzChunk + k] = -P( id, 3, i, j, k, index );
            // u[id * nzChunk + k ] =-P( id, 3, i, j, k, index );
            //      tmpMGReal[id * nzChunk + k + OFFS] =4.*scale*k;
            // cout<<" ( " << i <<" , "<< j <<" , "<< k<<" ) "<<tmpMGReal[k]<<endl;
            //     cout<<P( id, dir, i, j, k, index )<<endl;
        }
    }
}

#if ( PITTPACKACC )
#pragma acc routine gang
#endif
void MultiGrid::fillInArrayBack( ChunkedArray &P, int nChunk, int nzChunk, int i, int j, int index )
{
//     cout<<"after"<<endl;
#if ( PITTPACKACC )
#pragma acc loop gang
#endif
    for ( int id = 0; id < nChunk; id++ )
    {
#if ( PITTPACKACC )
#pragma acc loop vector
#endif
        for ( int k = 0; k < nzChunk; k++ )
        {
            P( id, 3, i, j, k, index ) = u[id * nzChunk + k]; //      cout<<k<<"   "<<tmpMG[i]<<endl;
            //    P( id, 3, i, j, k, index ) =5.*scale; //      cout<<k<<"   "<<tmpMG[i]<<endl;
            //        P(id, dir, i, j, k, 1 )  = 0; //      cout<<k<<"   "<<tmpMGReal[i]<<endl;
            //      cout<<tmpMGReal[k]<<endl;
        }
    }
    //      cout<<"-------------------"<<endl;
}

void MultiGrid::suggestSize( int a )
{
    int count = 0;

    for ( int i = 2; i < a / 2 + 1; i++ )
    {
        if ( a % i == 0 )
        {
            count++;
        }
    }

    // int *q=(int*)malloc(sizeof(int)*count);
    int *q = new int[count];

    count = 0;

    for ( int i = 2; i < a / 2 + 1; i++ )
    {
        if ( a % i == 0 )
        {
            q[count] = i;
            count++;
        }
    }

    cout << RED << " Suggestion " << a << " is divisible by: " << endl;
    for ( int i = 0; i < count; i++ )
    {
        cout << GREEN << q[i] << RESET << endl;
    }

    delete[] q;
}
