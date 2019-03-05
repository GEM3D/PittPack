#include "mathFunction.h"
#include <math.h>

#if ( OPENACC )
#include "openacc.h"
#endif

#if ( OPENACC )
#pragma acc routine seq
#endif
double sine( double x )
{
    return ( sin( x ) );
}

#if ( OPENACC )
#pragma acc routine seq
#endif
double cosine( double x )
{
    return ( cos( x ) );
}

#if ( OPENACC )
#pragma acc routine seq
#endif
double absolute( double x )
{
    return ( fabs( x ) );
}

#if ( OPENACC )
#pragma acc routine seq
#endif
double arctan( double x )
{
    return ( atan( x ) );
}

#if ( OPENACC )
#pragma acc routine seq
#endif
double squareRoot( double x )
{
    return ( sqrt( x ) );
}


#if ( OPENACC )
#pragma acc routine seq
#endif
int myLog2( int x ) 
{
    return ( log2( x ) );
}

