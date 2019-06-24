#ifndef _MATHFUNCTION_H_
#define _MATHFUNCTION_H_

#include "params.h"

#if ( PITTPACKACC )
#include <accel.h>
#else
#include <cmath>
#endif

#if ( PITTPACKACC )
#pragma acc routine seq
#endif
double sine( double x );

#if ( PITTPACKACC )
#pragma acc routine seq
#endif
double cosine( double x );

#if ( PITTPACKACC )
#pragma acc routine seq
#endif
double absolute( double x );

#if ( PITTPACKACC )
#pragma acc routine seq
#endif
double arctan( double x );

#if ( PITTPACKACC )
#pragma acc routine seq
#endif
double squareRoot( double x );

#if ( PITTPACKACC )
#pragma acc routine seq
#endif
int myLog2( int x );

#endif
