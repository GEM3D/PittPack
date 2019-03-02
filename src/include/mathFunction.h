#ifndef _MATHFUNCTION_H_
#define _MATHFUNCTION_H_

#if ( OPENACC )
//#include <accelmath.h>
#else
#include <cmath>
#endif

#if ( OPENACC )
#pragma acc routine seq
#endif
double sine( double x );

#if ( OPENACC )
#pragma acc routine seq
#endif
double cosine( double x );

#if ( OPENACC )
#pragma acc routine seq
#endif
double absolute( double x );

#if ( OPENACC )
#pragma acc routine seq
#endif
double arctan( double x );

#if ( OPENACC )
#pragma acc routine seq
#endif
double squareRoot(double x);

#if ( OPENACC )
#pragma acc routine seq
#endif
int mylog2(int x);

#endif
