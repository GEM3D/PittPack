#include "signalProc.hpp"
#include "chunkedArray.hpp"
#include "mathFunction.hpp"
#include <cmath>

void SignalProc::copyin()
{
#if ( PITTPACKACC )
#pragma acc enter data create( this [0:1] )
#pragma acc update device( this )
#endif
}

SignalProc::SignalProc()
{
#if ( 0 )
#pragma acc enter data create( this [0:1] )
#pragma acc update device( this )
#endif
}

SignalProc::~SignalProc()
{
#if ( PITTPACKACC )
#pragma acc exit data delete ( this )
#endif
}

#if ( PITTPACKACC )
#pragma acc routine vector
#endif
void SignalProc::postprocessSignalDCT10( ChunkedArray &P, const int size, const int i, const int direction )
{
    int idx = size % 2;

    // assign the conjugates
    if ( idx == 1 )
    {
        idx = 0;
    }
    else
    {
        idx = 1;
    }

    //    if ( direction == 0 )
    {
#if ( PITTPACKACC )
#pragma acc loop vector
#endif
        for ( int j = 0; j < size / 2 - idx; j++ )
        {
            P( 2 * ( i * size + j + size / 2 + 1 ) )     = P( 2 * ( i * size + size / 2 - j - idx ) );
            P( 2 * ( i * size + j + size / 2 + 1 ) + 1 ) = -P( 2 * ( i * size + size / 2 - j - idx ) + 1 );
        }
    }
    double theta;

// notice the multiplication by two here
#if ( PITTPACKACC )
#pragma acc loop vector
#endif
    for ( int j = 0; j < size; j++ )
    {
        theta = ( pi / 2. / size ) * j;

        P( 2 * ( i * size + j ) )
        = 2. * ( P( 2 * ( i * size + j ) ) * ( cosine( theta ) ) + P( 2 * ( i * size + j ) + 1 ) * sine( theta ) );
        P( 2 * ( i * size + j ) + 1 ) = 0.0;
    }
}

#if ( PITTPACKACC )
#pragma acc routine vector
#endif
void SignalProc::postprocessSignalDCT01( ChunkedArray &P, const int size, const int j, const int direction )
{
    // to prevent overwrite we use complex part and then swap real and complex
    // parts
    int upperBound = ( size - 1 ) / 2 + 1;
    //   int istart     = upperBound;

#if ( PITTPACKACC )
#pragma acc loop vector
#endif
    for ( int i = 0; i < upperBound; i++ )
    {
        P( 2 * ( j * size + 2 * i ) + 1 ) = P( 2 * ( j * size + i ) );
    }

    // int idx = size % 2;

#if ( PITTPACKACC )
#pragma acc loop vector
#endif
    for ( int i = 0; i < size / 2; i++ )
    {
        P( 2 * ( j * size + 2 * i + 1 ) + 1 ) = P( 2 * ( j * size + size - i - 1 ) );
    }
#if ( PITTPACKACC )
#pragma acc loop vector
#endif
    for ( int i = 0; i < size; i++ )
    {
        P( 2 * ( i + j * size ) )     = P( 2 * ( i + j * size ) + 1 );
        P( 2 * ( i + j * size ) + 1 ) = 0.0;
    }
}
#if ( PITTPACKACC )
#pragma acc routine vector
#endif
void SignalProc::preprocessSignalDCT01( ChunkedArray &P, const int size, const int i, const int direction )
{
    double theta0 = pi / 2. / size;

    double theta = 0.0;

    P( 2 * i * size + 0 ) = 0.5 * cosine( theta ) * P( 2 * i * size );
    P( 2 * i * size + 1 ) = 0.0;

#if ( PITTPACKACC )
#pragma acc loop vector
#endif
    for ( int j = 1; j < size / 2 + 1; j++ )
    {
        theta = theta0 * j;
        //
        // watch out the order, first assign the complex as I am performing in
        // place transform
        //
        P( 2 * ( i * size + j ) + 1 )
        = 0.5 * ( sine( theta ) * P( 2 * ( i * size + j ) ) - cosine( theta ) * P( 2 * ( i * size + size - j ) ) );
        P( 2 * ( i * size + j ) )
        = 0.5 * ( cosine( theta ) * P( 2 * ( i * size + j ) ) + sine( theta ) * P( 2 * ( i * size + size - j ) ) );
    }

#if ( PITTPACKACC )
#pragma acc loop vector
#endif
    for ( int j = size / 2 + 1; j < size; j++ )
    {
        P( 2 * ( i * size + j ) )     = P( 2 * ( i * size + size - j ) );
        P( 2 * ( i * size + j ) + 1 ) = -P( 2 * ( i * size + size - j ) + 1 );
    }

#if ( PITTPACKACC )
#pragma acc loop vector
#endif
    for ( int j = size / 2 + 1; j < size; j++ )
    {
        P( 2 * ( i * size + j ) )     = P( 2 * ( i * size + size - j ) );
        P( 2 * ( i * size + j ) + 1 ) = -P( 2 * ( i * size + size - j ) + 1 );
    }
}

#if ( PITTPACKACC )
//#pragma acc routine vector
#pragma acc routine vector
#endif
void SignalProc::preprocessSignalDCT10( ChunkedArray &P, const int size, const int i, const int direction )
{
    // cout<<RED<<P.ny<<" vs " << NXCHUNK  <<endl;

    int upperBound = ( size - 1 ) / 2 + 1;
    int istart     = upperBound;

#if ( PITTPACKACC )
#pragma acc loop vector
#endif
    for ( int j = 0; j < upperBound; j++ )
    {
        P( 2 * ( i * size + j ) + 1 ) = P( 2 * ( i * size + 2 * j ) );
    }

    int idx = size % 2;

#if ( PITTPACKACC )
#pragma acc loop vector
#endif
    for ( int j = istart; j < size; j++ )
    {
        P( 2 * ( i * size + j ) + 1 ) = P( 2 * i * size + 2 * ( size - 2 * ( j - istart ) - 1 - idx ) );
    }

    double tmp;
#if ( PITTPACKACC )
#pragma acc loop vector private( tmp )
#endif
    for ( int j = 0; j < size; j++ )
    {
        tmp                           = P( 2 * ( i * size + j ) + 1 );
        P( 2 * ( i * size + j ) )     = tmp;
        P( 2 * ( i * size + j ) + 1 ) = 0.0;
    }
}

// make sure the assignment from given real pointer P to copmex carray is consistent with the boundary  conditions
// this routine assumes that the signal is saved contigeously disregarding suxh that the signal fits in the first N/2 elements
#if ( PITTPACKACC )
#pragma acc routine vector
#endif
void SignalProc::preprocessSignalDST10( ChunkedArray &P, const int size, const int i, const int direction )
{
    swap( P, size, i, direction );
    preprocessSignalDCT10( P, size, i, direction );
}

#if ( PITTPACKACC )
#pragma acc routine vector
#endif
void SignalProc::postprocessSignalDST10( ChunkedArray &P, const int size, const int i, const int direction )
{
    postprocessSignalDCT10( P, size, i, direction );
    swap( P, size, i, direction );
}

#if ( PITTPACKACC )
#pragma acc routine vector
#endif
void SignalProc::preprocessSignalDST01( ChunkedArray &P, const int size, const int i, const int direction )
{
    swap( P, size, i, direction );
    preprocessSignalDCT01( P, size, i, direction );
}

#if ( PITTPACKACC )
#pragma acc routine vector
#endif
void SignalProc::postprocessSignalDST01( ChunkedArray &P, const int size, const int i, const int direction )
{
    postprocessSignalDCT01( P, size, i, direction );
    swap( P, size, i, direction );
}

// A DST-III  can be computed from a DCT-III or DCT-IV (see discrete cosine transform),
// respectively, by reversing the order of the inputs and flipping the sign of every other output,
// and vice versa for DST-II from DCT-II, ref : https://en.wikipedia.org/wiki/Discrete_sine_transform

#if ( PITTPACKACC )
#pragma acc routine vector
#endif
void SignalProc::swap( ChunkedArray &P, const int size, const int i, const int direction )
{
#if ( PITTPACKACC )
#pragma acc loop vector
#endif
    for ( int j = 0; j < size; j++ )
    {
        P( 2 * ( i * size + j ) + 1 ) = P( 2 * ( i * size + size - j - 1 ) );

        if ( j % 2 == 1 )
        {
            P( 2 * ( i * size + j ) + 1 ) = -P( 2 * ( i * size + j ) + 1 );
        }
    }

#if ( PITTPACKACC )
#pragma acc loop vector
#endif
    for ( int j = 0; j < size; j++ )
    {
        P( 2 * ( i * size + j ) )     = P( 2 * ( i * size + j ) + 1 );
        P( 2 * ( i * size + j ) + 1 ) = 0.0;
    }
}
