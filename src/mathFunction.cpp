#include "mathFunction.hpp"
#include <math.h>

#if ( PITTPACKACC )
#include "openacc.h"
#endif

#if ( PITTPACKACC )
#pragma acc routine seq
#endif
double sine( double x ) { return ( sin( x ) ); }

#if ( PITTPACKACC )
#pragma acc routine seq
#endif
double cosine( double x ) { return ( cos( x ) ); }

#if ( PITTPACKACC )
#pragma acc routine seq
#endif
double absolute( double x ) { return ( fabs( x ) ); }

#if ( PITTPACKACC )
#pragma acc routine seq
#endif
double arctan( double x ) { return ( atan( x ) ); }

#if ( PITTPACKACC )
#pragma acc routine seq
#endif
double squareRoot( double x ) { return ( sqrt( x ) ); }

#if ( PITTPACKACC )
#pragma acc routine seq
#endif
int myLog2( int x ) { return ( log2( x ) ); }
