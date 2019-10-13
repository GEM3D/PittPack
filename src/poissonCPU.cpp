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

#if 0
void PoissonCPU::pittPack()
{
#if ( DEBUG0 )
    ofstream myfile;

    std::string filename = "data";
    filename.append( to_string( myRank ) );
    //  ofstream myfile;
    myfile.open( filename );
#endif

    //    MPI_Barrier( MPI_COMM_WORLD );
    // double t1 = MPI_Wtime();

    double t1 = 0.0;
    double t2 = 0.0;

    //  struct timeval start_time, stop_time, elapsed_time;
    //   gettimeofday( &start_time, NULL );
    int result;

    double err    = 0.0;
    double t1_com = 0.0;
    double t2_com = 0.0;
    double deT    = 0.0;

    for ( int num = 0; num < NITER; num++ )
    // for(int num=0;num<100;num++)
    {
        if ( INITANALYTIC )
        {
            initializeTrigonometric();
        }
        //   if ( bc[0] != 'P' && bc[2] != 'P' )
#if ( SOLVE != 0 )
        modifyRhsDirichlet();
#endif

        t1 = MPI_Wtime();

#if ( POSS )

        if ( PROFILE_COMM )
        {
            t1_com = MPI_Wtime();
            changeOwnershipPairwiseExchangeZX();
            t2_com = MPI_Wtime() - t1_com;
        }
        else
        {
            changeOwnershipPairwiseExchangeZX();
        }

        //    M.rearrange( 0, 0, 1 );
        //   M.rearrange( 0, 2 );

#if ( DEBUG0 )
        myfile << "     Z to X rotation" << endl;
        printX( myfile );
#endif

        // perform FFT in x direction
        // remember, need to rearrange and restore each time

        // step 2) change location of the array such that FFT can be performed on a contegeous array

#if ( COMM_PATTERN == 2 )
        changeLocationXOverlap();
#else
        changeLocationX();
#endif

#if ( DEBUG0 )
        myfile << "     change Loc in X-- where FFT should be envoked" << endl;
        printX( myfile );
#endif
        // step 3) perform  FFT
        preprocessSignalAccordingly( 0, 0 );

#if ( DEBUG0 )
        myfile << "     preprocess Signal before FFTX" << endl;
        // M.printX( myfile );
        printX( myfile );

#endif

        performTransformXdir();

#if ( FFTX )
#if ( DEBUG0 )
        myfile << "     FFTX" << endl;
        printX( myfile );

#endif

        postprocessSignalAccordingly( 0, 0 );

#if ( DEBUG0 )
        myfile << "     postprocess Signal after FFTX" << endl;
        // M.printX( myfile );
        printX( myfile );

#endif

        // step 4) restore the array to original status before FFT

        restoreLocationX();
#if ( DEBUG0 )
        myfile << "     Restore Loc in X" << endl;
        printX( myfile );
#endif
#endif

        //       cout << "FFTX is done" << endl;
#if ( FFTY )
        // step 5) pencils with n(1,0,0) is converted to pencil with n(0,1,0)

        if ( PROFILE_COMM )
        {
            t1_com = MPI_Wtime();
            changeOwnershipPairwiseExchangeXY();
            t2_com += MPI_Wtime() - t1_com;
        }
        else
        {
            changeOwnershipPairwiseExchangeXY();
        }

#if ( DEBUG0 )
        myfile << "     X to Y rotation" << endl;
        printY( myfile );
#endif
        //    M.rearrange( 0, 1 );
        // step 6) swaps X and Y coordinates, now X is Y and nx is ny

        rearrangeX2Y();

#if ( DEBUG0 )
        myfile << "     Rearrange each block from X to Y direction" << endl;

        // M.printX( myfile );
        printY( myfile );
#endif

        // step 7) change location of the array such that FFT can be performed on a contegeous array in the transverse direction

        changeLocationY();

#if ( DEBUG0 )
        myfile << "      change Loc in Y-- where FFT should be envoked in Y-direction" << endl;
        // M.printX( myfile );
        printY( myfile );

// step 8) perform  FFT transform can be performed
#endif
        preprocessSignalAccordingly( 1, 1 );

#if ( DEBUG0 )
        myfile << "     preprocess Signal before FFTY" << endl;
        // M.printX( myfile );
        printY( myfile );

// step 8) perform  FFT transform can be performed
#endif

        performTransformYdir();
        // M.printX( myfile );

#if ( DEBUG0 )
        myfile << "     FFTY" << endl;
        printY( myfile );

#endif

        postprocessSignalAccordingly( 1, 1 );
#if ( DEBUG0 )
        myfile << "     postprocess Signal before FFTY" << endl;
        // M.printX( myfile );
        printY( myfile );
#endif

        // step 9) restore the array to original status before FFT

        restoreLocationY();
#if ( DEBUG0 )
        myfile << "      Restore Loc in Y" << endl;
        // M.printX( myfile );
        printY( myfile );
#endif
#endif
        // cout << "FFTY is done" << endl;
        // step 10) pencils with n(0,1,0) is converted to pencil with n(0,0,1)

#if ( SOLVE )

        if ( PROFILE_COMM )
        {
            t1_com = MPI_Wtime();
            //   M.changeOwnershipPairwiseExchangeYZ();
            changeOwnershipPairwiseExchangeZX();
            t2_com += MPI_Wtime() - t1_com;
        }
        else
        {
            changeOwnershipPairwiseExchangeZX();
        }

#if ( DEBUG0 )
        myfile << "      Rotate Y to Z " << endl;
        // M.printX( myfile );
        printY( myfile );
        myfile << "      EigenValue" << endl;
        eigenVal( myfile );
#endif
        // step 11) Customized multiBlock Thomas and periodic Thomas (with Sherman-Morrison modification)

        //
        if ( SOLUTIONMETHOD == 0 || SOLUTIONMETHOD == 2 )
        {
            /*
                        if ( solveThm( 0 ) != SUCCESS )
                        {

                            solveThm( 1 );
                            cout << "Exit Code: " << THOMAS_FAIL << endl;
                            cout << BLUE << PittPackGetErrorEnum( THOMAS_FAIL ) << RESET << endl;
                            exit( 1 );


                        }

            */

            solveThmBatch( 0 );

            if ( bc[0] == 'P' || bc[2] == 'P' )
            {
                solveThmBatch( 1 );
            }

            //         solveThm(0);

            //         solveCRP(0);
        }
        else if ( SOLUTIONMETHOD == 1 )
        {
            // "0" solves for the real part and "1" solves for imaginary part
            // solveMG();
            //   solveMGC( );
            //        solveCRP(0);
        }
        // M.printX( myfile );
        //

#if ( DEBUG0 )
        myfile << "      Thomasing" << endl;
        printY( myfile );

#endif
        // step 12) pencils with n(0,0,1) is converted to pencil with n(0,1,0)
        // M.changeOwnershipPairwiseExchangeYZ();
        //

        if ( PROFILE_COMM )
        {
            t1_com = MPI_Wtime();
            changeOwnershipPairwiseExchangeZX();
            t2_com += MPI_Wtime() - t1_com;
        }
        else
        {
            changeOwnershipPairwiseExchangeZX();
        }

#if ( DEBUG0 )
        myfile << "      Rotate Z to Y" << endl;
        // M.printX( myfile );
        printY( myfile );
#endif
#endif

#if ( IFFTY )

        // step 13) prepre for contigeuous FFT

        changeLocationY();
#if ( DEBUG0 )
        // M.printX( myfile );
        myfile << "     change Loc in Y to perform IFFTY" << endl;
        printY( myfile );
#endif
        // step 14) perform IFFT

        preprocessSignalAccordinglyReverse( 1, 1 );

#if ( DEBUG0 )
        myfile << "     preprocess Signal before IFFTY" << endl;
        // M.printX( myfile );
        printY( myfile );

#endif

        performInverseTransformYdir();
#if ( DEBUG0 )
        myfile << "     perform IFFTY" << endl;
        printY( myfile );
#endif

        postprocessSignalAccordinglyReverse( 1, 1 );
// step 15) restore the array to original status before IFFT
#if ( DEBUG0 )
        myfile << "     postprocess Signal before IFFTY" << endl;

        printY( myfile );
#endif

        restoreLocationY();

#if ( DEBUG0 )
        myfile << "      Restore Loc in Y" << endl;
        // M.printX( myfile );
        printX( myfile );
#endif

        // step 16) swaps X and Y coordinates, now X is Y and nx is ny
        // need to redefine the function ??????????????????????????????????????????????????????????????????????????
        // now again switching to x-direction arrangement so printX is fine
        rearrangeX2YInverse();

#if ( DEBUG0 )
        myfile << "     Rearrange each block from Y to X direction" << endl;
        printX( myfile );
#endif
        // step 17) pencils with n(0,1,0) is converted to pencil with n(1,0,0)

        if ( PROFILE_COMM )
        {
            t1_com = MPI_Wtime();
            changeOwnershipPairwiseExchangeXY();
            t2_com += MPI_Wtime() - t1_com;
        }
        else
        {
            changeOwnershipPairwiseExchangeXY();
        }
#if ( DEBUG0 )
        myfile << "      Rotate Y to X" << endl;
        printX( myfile );

#endif
#endif

#if ( IFFTX )
        // step 18) change location of the array such that IFFT can be performed on a contegeous array
        changeLocationX();

#if ( DEBUG0 )
        myfile << "     change Loc in X to perform IFFTX" << endl;

        printX( myfile );

#endif
        // step 19) perform  IFFT

        preprocessSignalAccordinglyReverse( 0, 0 );
#if ( DEBUG0 )
        myfile << "     preprocess Signal before IFFTX" << endl;
        printX( myfile );
#endif

        performInverseTransformXdir();

#if ( DEBUG0 )
        myfile << "     perform IFFTX" << endl;
        printX( myfile );
#endif

        postprocessSignalAccordinglyReverse( 0, 0 );
// step 20) restore the array to original status before FFT
#if ( DEBUG0 )
        myfile << "     postprocess Signal before IFFTY" << endl;
        printX( myfile );
#endif

        restoreLocationX();
#if ( DEBUG0 )
        myfile << "     Final result" << endl;
        printX( myfile );
#endif
        rescale();

        t2 = MPI_Wtime();

        deT += ( t2 - t1 );

#if ( DEBUG0 )
        myfile << "     Final result Rescaled" << endl;
        printX( myfile );

#endif

#endif

        // return back to the original set-up
        if ( PROFILE_COMM )
        {
            t1_com = MPI_Wtime();
            changeOwnershipPairwiseExchangeZX();
            t2_com += MPI_Wtime() - t1_com;
        }
        else
        {
            changeOwnershipPairwiseExchangeZX();
        }

        t2 = MPI_Wtime();
        deT += ( t2 - t1 );

        if ( INCLUDE_ERROE_CAL_IN_TIMING == 1 )
        {
            err = getError();
        }

#endif
    }
    cout << " Running .......  " << endl;

    if ( INCLUDE_ERROE_CAL_IN_TIMING == 1 )
    {
        MPI_Allreduce( &err, &finalErr, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
        finalErr = finalErr / p0 / p0;

        if ( myRank == 0 )
        {
            // cout << YELLOW << "Error  (" << myRank << ") =" << err << " " << finalErr << RESET << endl;
            cout << YELLOW << "Error Per processor"
                 << " " << finalErr << RESET << endl;
        }
    }

#if ( DEBUG0 )
    myfile.close();
#endif

    MPI_Reduce( &deT, &runTime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD );

    if ( myRank == 0 )
    {
        runTime = runTime / (double)NITER;
        if(RUNINFO)
        {
        runInfo();
        }
        if ( PROFILE_COMM )
        {
            printf( "change ownership time =%lf percent of solution and take %lf seconds \n", t2_com / deT * 100., t2_com );
        }
    }

    if ( ZERO_MEAN )
    {
        subtractMeanValue();
    }
}

#endif
/* absolute value version
void PoissonCPU::testDST10()
{
    int mysize = 4;

    fftw_plan pl;

    double *inEx = new double[2*(mysize)];
    double *in = new double[2*(mysize)];
    double *out = new double[2*mysize];
//    fftw_complex *outC = (fftw_complex *)malloc( sizeof( fftw_complex ) * 2 * ( mysize + 1 ) );
//    fftw_complex *inC = (fftw_complex *)malloc( sizeof( fftw_complex ) * 2 * ( mysize + 1 ) );


// to obtain DST2
    for ( int i = 0; i < mysize; i++ )
    {
        // in[i]=exp(i);
    //    inC[i][0] = ( i ) * i;
    //    inC[i][1] = 0.0;
        in[i] = (( i ) * (i));
     //      in[i] = rand() % 200000;
    }

    int size = mysize;

    pl = fftw_plan_r2r_1d( size, in, out, FFTW_RODFT10, FFTW_ESTIMATE );

    fftw_execute( pl );

    cout << "original" << endl;

    for ( int i = 0; i < 2*mysize; i++ )
    {
        cout << in[i] << endl;
    }

    cout << "real transform" << endl;

    for ( int i = 0; i < mysize; i++ )
    {
        cout << out[i] << endl;
    }


    fftw_complex *outC = (fftw_complex *)malloc( sizeof( fftw_complex ) * (2* mysize  ) );
    fftw_complex *inC = (fftw_complex *)malloc( sizeof( fftw_complex ) * 2 * ( mysize  ) );

   for ( int i = 0; i < mysize; i++ )
    {
        // in[i]=exp(i);
    //    inC[i][0] = ( i ) * i;
    //    inC[i][1] = 0.0;
       in[mysize+i] =-in[mysize-i-1];
    }
  cout << "complex transform" << endl;
   for ( int i = 0; i < mysize; i++ )
    {
        // in[i]=exp(i);
        inC[i][0] =- in[i];
        inC[i][1]=0.0;

    //    inC[i][1] = 0.0;
     //   in[mysize+i] =in[mysize-i-1];
        cout<<inC[i][0]<<"+ "<<inC[i][1]<<"i"<<endl;
    }



    pl = fftw_plan_dft_r2c_1d( 2*size, in, outC, FFTW_ESTIMATE );
//    pl = fftw_plan_dft_1d( 2*size, inC, outC,FFTW_FORWARD, FFTW_ESTIMATE );

    fftw_execute( pl );

  cout << "complex transform" << endl;

    for ( int i = 0; i <mysize+1; i++ )
    {
        cout << outC[i][0] <<"+ "<<outC[i][1]<<"i "<< endl;
    }




}
*/
/*
void PoissonCPU::testDST01()
{
    int mysize =7 ;

    fftw_plan pl;

    double *inEx = new double[2*(mysize)];
    double *in = new double[2*(mysize)];
    double *out = new double[2*mysize];
//    fftw_complex *outC = (fftw_complex *)malloc( sizeof( fftw_complex ) * 2 * ( mysize + 1 ) );
//    fftw_complex *inC = (fftw_complex *)malloc( sizeof( fftw_complex ) * 2 * ( mysize + 1 ) );


// to obtain DST2
    for ( int i = 0; i < mysize; i++ )
    {
        // in[i]=exp(i);
    //    inC[i][0] = ( i ) * i;
    //    inC[i][1] = 0.0;
        in[i] = (( i+1 ) * (i+1));
    }

    int size = mysize;

    pl = fftw_plan_r2r_1d( size, in, out, FFTW_RODFT01, FFTW_ESTIMATE );

    fftw_execute( pl );

    cout << "original" << endl;

    for ( int i = 0; i < 2*mysize; i++ )
    {
        cout << in[i] << endl;
    }

    cout << "real transform" << endl;

    for ( int i = 0; i < mysize; i++ )
    {
        cout << out[i] << endl;
    }


    fftw_complex *outC = (fftw_complex *)malloc( sizeof( fftw_complex ) * (2* mysize  ) );
    fftw_complex *inC = (fftw_complex *)malloc( sizeof( fftw_complex ) * 2 * ( mysize  ) );

   for ( int i = 1; i < mysize; i++ )
    {
        // in[i]=exp(i);
    //    inC[i][0] = ( i ) * i;
    //    inC[i][1] = 0.0;
       in[mysize+i-1] =in[mysize-i-1];
    }
  cout << "complex transform" << endl;
//  inC[0][0]=0.0;
   for ( int i = 0; i < 2*mysize; i++ )
    {
        // in[i]=exp(i);
        inC[i][0] = in[i];
        inC[i][1]=0.0;

    //    inC[i][1] = 0.0;
     //   in[mysize+i] =in[mysize-i-1];
        cout<<inC[i][0]<<"+ "<<inC[i][1]<<"i"<<endl;
    }



//      pl = fftw_plan_dft_r2c_1d( 2*size, in, outC, FFTW_ESTIMATE );
//    pl = fftw_plan_dft_1d( 2*size, inC, outC,FFTW_FORWARD, FFTW_ESTIMATE );

    fftw_execute( pl );

  cout << "complex transform" << endl;

    for ( int i = 0; i < 2*mysize; i++ )
    {
        cout << outC[i][0] <<"+ "<<outC[i][1]<<"i "<< endl;
    }




}
*/
#if 0 
// It is possible to calculate the DST01 from DCT01 and yes DCT01, there are no Typos :-)

void PoissonCPU::testDST01()
{
    int mysize = 7;

    fftw_plan pl;

    double *inEx = new double[2 * ( mysize )];
    double *in   = new double[2 * ( mysize )];
    double *out  = new double[2 * mysize];

    // to obtain DST2
    for ( int i = 0; i < mysize; i++ )
    {
        in[i] = ( ( i + 1 ) * ( i + 1 ) );
    }

    int size = mysize;

    pl = fftw_plan_r2r_1d( size, in, out, FFTW_RODFT01, FFTW_ESTIMATE );

    fftw_execute( pl );

    cout << "original" << endl;

    for ( int i = 0; i < mysize; i++ )
    {
        cout << in[i] << endl;
    }

    cout << "real transform" << endl;

    for ( int i = 0; i < mysize; i++ )
    {
        cout << out[i] << endl;
    }

    double *inC  = new double[2 * ( mysize )];
    double *outC = new double[2 * mysize];

    for ( int i = 0; i < mysize; i++ )
    {
        inC[i] = in[mysize - i - 1];

        if ( i % 2 == 1 )
        {
            inC[i] = -inC[i];
        }
    }

    pl = fftw_plan_r2r_1d( size, inC, out, FFTW_REDFT01, FFTW_ESTIMATE );

    cout << "complex transform" << endl;

    for ( int i = 0; i < mysize; i++ )
    {
        cout << inC[i] << endl;
    }

    fftw_execute( pl );
    /*
       for ( int i = 0; i < mysize; i++ )
          {
             out[i] =out[mysize-i-1];

           if(i%2==1)
           {
           out[i]=-out[i];
            }

        }
    */
 fftw_destroy_plan(pl);
    fftw_cleanup();


    cout << "complex transform" << endl;

    for ( int i = 0; i < mysize; i++ )
    {
        //       cout << outC[i][0] <<"+ "<<outC[i][1]<<"i "<< endl;
        cout << out[i] << endl;
    }
}

void PoissonCPU::testDST10()
{
    int mysize = 7;

    fftw_plan pl;

    double *inEx = new double[2 * ( mysize )];
    double *in   = new double[2 * ( mysize )];
    double *out  = new double[2 * mysize];

    // to obtain DST2
    for ( int i = 0; i < mysize; i++ )
    {
        in[i] = ( ( i + 1 ) * ( i + 1 ) );
    }

    int size = mysize;

    pl = fftw_plan_r2r_1d( size, in, out, FFTW_RODFT10, FFTW_ESTIMATE );

    fftw_execute( pl );

    cout << "original" << endl;

    for ( int i = 0; i < mysize; i++ )
    {
        cout << in[i] << endl;
    }

    cout << "real transform" << endl;

    for ( int i = 0; i < mysize; i++ )
    {
        cout << out[i] << endl;
    }

    double *inC  = new double[2 * ( mysize )];
    double *outC = new double[2 * mysize];

    for ( int i = 0; i < mysize; i++ )
    {
        inC[i] = in[mysize - i - 1];

        if ( i % 2 == 1 )
        {
            inC[i] = -inC[i];
        }
    }

    pl = fftw_plan_r2r_1d( size, inC, out, FFTW_REDFT10, FFTW_ESTIMATE );

    cout << "complex transform" << endl;

    for ( int i = 0; i < mysize; i++ )
    {
        cout << inC[i] << endl;
    }

    fftw_execute( pl );

    /*  rearrnage needs a coontainer
       for ( int i = 0; i < mysize; i++ )
          {
             out[i] =out[mysize-i-1];

           if(i%2==1)
           {
           out[i]=-out[i];
            }

        }
    */

    cout << "complex transform" << endl;

    for ( int i = 0; i < mysize; i++ )
    {
        //       cout << outC[i][0] <<"+ "<<outC[i][1]<<"i "<< endl;
        cout << out[i] << endl;
    }
 fftw_destroy_plan(pl);
    fftw_cleanup();


}

void checkCudaSupport()
{
    printf( "Compile time check:\n" );
#if defined( MPIX_CUDA_AWARE_SUPPORT ) && MPIX_CUDA_AWARE_SUPPORT
    printf( "This MPI library has CUDA-aware support.\n", MPIX_CUDA_AWARE_SUPPORT );
#elif defined( MPIX_CUDA_AWARE_SUPPORT ) && !MPIX_CUDA_AWARE_SUPPORT
    printf( "This MPI library does not have CUDA-aware support.\n" );
#else
    printf( "This MPI library cannot determine if there is CUDA-aware support.\n" );
#endif /* MPIX_CUDA_AWARE_SUPPORT */

    printf( "Run time check:\n" );
#if defined( MPIX_CUDA_AWARE_SUPPORT )
    if ( 1 == MPIX_Query_cuda_support() )
    {
        printf( "This MPI library has CUDA-aware support.\n" );
    }
    else
    {
        printf( "This MPI library does not have CUDA-aware support.\n" );
    }
#else  /* !defined(MPIX_CUDA_AWARE_SUPPORT) */
    printf( "This MPI library cannot determine if there is CUDA-aware support.\n" );
#endif /* MPIX_CUDA_AWARE_SUPPORT */

    //#pragma acc update self(a,b)

    // cout<<b  <<endl;
}



#endif
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

