#include "definitions.h"
#include "params.h"
#include "pencilDcmp.hpp"

void PoissonCPU::writeYLine( int j, fftw_complex *outC )
{
#if(!PITTPACKACC)

    for ( int i = 0; i < nChunk * nyChunk; i++ )
    {
        PencilDcmp::P( 2 * nyChunk * nChunk * j + 2 * i )     = outC[i][0];
        PencilDcmp::P( 2 * nyChunk * nChunk * j + 2 * i + 1 ) = outC[i][1];
#if ( DEBUG0 )
        cout << outC[i][0] << " " << outC[i][1] << '\t';
#endif
    }
#if ( DEBUG0 )
    cout << endl;
#endif
#endif
}

void PoissonCPU::readYLine( int j, fftw_complex *out )
{
#if(!PITTPACKACC)
    for ( int i = 0; i < nChunk * nyChunk; i++ )
    {
        out[i][0] = PencilDcmp::P( 2 * ( nChunk * nyChunk ) * j + 2 * i );
        out[i][1] = PencilDcmp::P( 2 * ( nChunk * nyChunk ) * j + 2 * i + 1 );
#if ( DEBUG0 )
        cout << out[i][0] << " " << out[i][1] << '\t';
#endif
    }
#if ( DEBUG0 )
    cout << endl;
#endif
#endif
}

void PoissonCPU::performTransformYdir()
{
//cout<<" cpu version "<<endl;
#if(!PITTPACKACC)
    fftw_plan pl;

    fftw_complex *in;
    in = (fftw_complex *)fftw_malloc( nChunk * nyChunk * sizeof( fftw_complex ) );

    fftw_complex *out;
    out = (fftw_complex *)fftw_malloc( nChunk * nyChunk * sizeof( fftw_complex ) );

    for ( int j = 0; j < nxChunk * nzChunk; j++ )
    {
        readYLine( j, in );

        pl = fftw_plan_dft_1d( nyChunk * nChunk, in, out, FFTW_FORWARD, FFTW_ESTIMATE );

        fftw_execute( pl );
        
        fftw_destroy_plan(pl);
        fftw_cleanup();
        writeYLine( j, out );
    }

    fftw_free( out );
    fftw_free( in );
#endif
//    fftw_destroy_plan(pl);
//    fftw_cleanup();

}

void PoissonCPU::performInverseTransformYdir()
{
#if(!PITTPACKACC)
    fftw_plan     pl;
    fftw_complex *in;
    in = (fftw_complex *)fftw_malloc( nChunk * nyChunk * sizeof( fftw_complex ) );

    fftw_complex *out;
    out = (fftw_complex *)fftw_malloc( nChunk * nyChunk * sizeof( fftw_complex ) );

    for ( int j = 0; j < nxChunk * nzChunk; j++ )
    {
        readYLine( j, in );

        pl = fftw_plan_dft_1d( nyChunk * nChunk, in, out, FFTW_BACKWARD, FFTW_ESTIMATE );

        fftw_execute( pl );
       fftw_destroy_plan(pl);
       fftw_cleanup();

       writeYLine( j, out );
    }

    
    fftw_free( out );
    fftw_free( in );
#endif

}

void PoissonCPU::readXLine( int j, fftw_complex *out )
{
#if(!PITTPACKACC)
    for ( int i = 0; i < nChunk * nxChunk; i++ )
    {
        out[i][0] = PencilDcmp::P( 2 * ( nChunk * nxChunk ) * j + 2 * i );
        out[i][1] = PencilDcmp::P( 2 * ( nChunk * nxChunk ) * j + 2 * i + 1 );
#if ( DEBUG0 )
        if ( This == myRank )
        {
            cout << GREEN << " " << out[i][0] << " " << out[i][1] << '\t';
        }
#endif
    }

#if ( DEBUG0 )
    cout << RESET << endl;
#endif
#endif
}

void PoissonCPU::writeXLine( int j, fftw_complex *outC )
{
#if(!PITTPACKACC)
    for ( int i = 0; i < nChunk * nxChunk; i++ )
    {
        PencilDcmp::P( 2 * nxChunk * nChunk * j + 2 * i )     = outC[i][0];
        PencilDcmp::P( 2 * nxChunk * nChunk * j + 2 * i + 1 ) = outC[i][1];
//     cout << " ???????????????????????" << out[i] << endl;
#if ( DEBUG0 )
        if ( This == myRank )
        {
            cout << BLUE << outC[i][0] << " " << outC[i][1] << '\t';
        }
#endif
    }
#if ( DEBUG0 )
    if ( This == myRank )
    {
        cout << RESET << endl;
    }
#endif
#endif
}

void PoissonCPU::performInverseTransformXdir()
{
#if(!PITTPACKACC)
    fftw_plan     pl;
    fftw_complex *in;
    in = (fftw_complex *)fftw_malloc( nChunk * nxChunk * sizeof( fftw_complex ) );

    fftw_complex *out;
    out = (fftw_complex *)fftw_malloc( nChunk * nxChunk * sizeof( fftw_complex ) );

    for ( int j = 0; j < nyChunk * nzChunk; j++ )
    {
        readXLine( j, in );

        pl = fftw_plan_dft_1d( nxChunk * nChunk, in, out, FFTW_BACKWARD, FFTW_ESTIMATE );

        fftw_execute( pl );
   fftw_destroy_plan(pl);
    fftw_cleanup();


        writeXLine( j, out );
    }

    fftw_free( out );
    fftw_free( in );
 
#endif
}

void PoissonCPU::performTransformXdir()
{
//cout <<" here "<<endl;
#if(!PITTPACKACC)
    fftw_plan pl;

    fftw_complex *in;
    in = (fftw_complex *)fftw_malloc( nChunk * nxChunk * sizeof( fftw_complex ) );

    if ( in == NULL )
    {
        cout << "error oin allocation " << endl;
    }

    fftw_complex *out;
    out = (fftw_complex *)fftw_malloc( nChunk * nxChunk * sizeof( fftw_complex ) );

    if ( out == NULL )
    {
        cout << "error oin allocation " << endl;
    }

    for ( int j = 0; j < nyChunk * nzChunk; j++ )
    {
#if ( DEBUG0 )
        if ( This == myRank )
        {
            cout << " in & myrank " << myRank << " " << j << endl;
        }
#endif

        readXLine( j, in );

        pl = fftw_plan_dft_1d( nxChunk * nChunk, in, out, FFTW_FORWARD, FFTW_ESTIMATE );
#if ( DEBUG0 )
        if ( This == myRank )
        {
            cout << " out " << j << endl;
        }
#endif

        fftw_execute( pl );
    fftw_destroy_plan(pl);
    fftw_cleanup();


        writeXLine( j, out );
    }

    //
    // cout<<RED<<nChunk << " "<<nxChunk<<RESET<<endl;
    fftw_free( out );
    fftw_free( in );

#endif
//  cout<< " poisson is calling "<<endl;
}


/*
    // always check for compatibility of the Class on GPU
    if ( my_rank == 0 )
    {
        if ( __has_trivial_copy( PoissonGPU ) == true )
        {
            std::cout << "Congrat! myClass will work with OpenACC copy operations " << std::endl;
        }
        else
        {
            std::cout << "OOPS! myClass will NOT work with OpenACC copy operations " << std::endl;
            return 1;
        }
    }
*/

