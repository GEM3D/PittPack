#include "definitions.h"
#include "pencilDcmp.hpp"

#if ( PITTPACKACC )
void PoissonGPU::performTransformXdir() /*!< Is called on Host and Runs on GPU*/
{
    cufftHandle plan = NULL;
    double *    ptr  = P.P;

#pragma acc host_data use_device( ptr )
    {
        cufftPlan1d( &plan, nx, CUFFT_Z2Z, nyChunk * nzChunk );
        if ( CUFFT_STREAMS )
        {
            void *stream = acc_get_cuda_stream( acc_async_sync );
            cufftSetStream( plan, (cudaStream_t)stream );
        }
        if ( CUFFT_SUCCESS != cufftExecZ2Z( plan, (cufftDoubleComplex *)( ptr ), (cufftDoubleComplex *)( ptr ), CUFFT_FORWARD ) )
        {
            cout << " Exit Code : " << PittPackGetErrorEnum( CUFFT_FAIL_X ) << endl;
            exit( CUFFT_FAIL_X );
        }
        cufftDestroy( plan );
    }
}

void PoissonGPU::performInverseTransformXdir() /*!< Called on Host and Ran on GPU*/
{
    cufftHandle plan = NULL;
    double *    ptr  = P.P;

#pragma acc host_data use_device( ptr )
    {
        cufftPlan1d( &plan, nx, CUFFT_Z2Z, nyChunk * nzChunk );
        if ( CUFFT_STREAMS )
        {
            void *stream = acc_get_cuda_stream( acc_async_sync );
            cufftSetStream( plan, (cudaStream_t)stream );
        }
        if ( CUFFT_SUCCESS != cufftExecZ2Z( plan, (cufftDoubleComplex *)( ptr ), (cufftDoubleComplex *)( ptr ), CUFFT_INVERSE ) )
        {
            cout << " Exit Code : " << PittPackGetErrorEnum( CUFFT_FAIL_INV_X ) << endl;
            exit( CUFFT_FAIL_INV_X );
        }
        cufftDestroy( plan );
    }
}
void PoissonGPU::performTransformYdir() /*!< Called on Host and Ran on GPU*/
{
    cufftHandle plan = NULL;
    double *    ptr  = P.P;

#pragma acc host_data use_device( ptr )
    {
        cufftPlan1d( &plan, ny, CUFFT_Z2Z, nxChunk * nzChunk );
        if ( CUFFT_STREAMS )
        {
            void *stream = acc_get_cuda_stream( acc_async_sync );
            cufftSetStream( plan, (cudaStream_t)stream );
        }
        if ( CUFFT_SUCCESS != cufftExecZ2Z( plan, (cufftDoubleComplex *)( ptr ), (cufftDoubleComplex *)( ptr ), CUFFT_FORWARD ) )
        {
            cout << " Exit Code : " << PittPackGetErrorEnum( CUFFT_FAIL_Y ) << endl;
            exit( CUFFT_FAIL_Y );
        }
        cufftDestroy( plan );
    }
}

void PoissonGPU::performInverseTransformYdir() /*!< Called on Host and Ran on GPU*/
{
    cufftHandle plan = NULL;
    double *    ptr  = P.P;

#pragma acc host_data use_device( ptr )
    {
        cufftPlan1d( &plan, ny, CUFFT_Z2Z, nxChunk * nzChunk );
        if ( CUFFT_STREAMS )
        {
            void *stream = acc_get_cuda_stream( acc_async_sync );
            cufftSetStream( plan, (cudaStream_t)stream );
        }
        if ( CUFFT_SUCCESS != cufftExecZ2Z( plan, (cufftDoubleComplex *)( ptr ), (cufftDoubleComplex *)( ptr ), CUFFT_INVERSE ) )
        {
            cout << " Exit Code : " << PittPackGetErrorEnum( CUFFT_FAIL_INV_Y ) << endl;
            exit( CUFFT_FAIL_INV_Y );
        }
        cufftDestroy( plan );
    }
}

void PoissonGPU::triDiagCusparse( double *dl, double *ds, double *du, double *rhs )
{
    // cout<< " size of nz" <<nz<<endl;

    cusparseHandle_t handle = NULL;

#pragma acc host_data use_device( rhs, du, dl, ds )
    {
        cusparseCreate( &handle );
        if ( PIVOT == 0 )
        {
            cusparseDgtsv_nopivot( handle, nz, 1, dl, ds, du, rhs, nz );
        }
        else
        {
            cusparseDgtsv( handle, nz, 1, dl, ds, du, rhs, nz );
        }
    }

    cusparseDestroy( handle );
}

// const std::string currentDateTime() ;

void PoissonGPU::pittPack() /*!<called on CPU runs on GPU */
{
    double t1 = 0.0;
    double t2 = 0.0;

#if ( DEBUG2 )
    ofstream myfile;

    std::string filename = "data";
    filename.append( to_string( myRank ) );
    //  ofstream myfile;
    myfile.open( filename );

    myfile << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx" << endl;
    myfile << "     Start" << endl;

    printX( myfile );

#endif

    //   initializeAndBind();

    long int initMem = acc_get_memory() / 1e9;

    if ( MONITOR_MEM )
    {
        cout << "amount of available memory " << initMem << endl;
    }
    double err = 0.0;
    finalErr   = 0.0;

    // double eig;

    double t1_com = 0.0;
    double t2_com = 0.0;
    double deT    = 0.0;

    //    P.moveHostToDevice();

#if ( 1 )

    int trsps_gang0 = MIN( 1024, iaxSize );
    int trsps_gang1 = MIN( 1024, iaySize );

    int nSig0 = nxChunk;
    int nSig1 = nyChunk;
    initMem   = acc_get_free_memory() / 1e9;

    if ( MONITOR_MEM )
    {
        cout << "amount of free memory before data region " << acc_get_free_memory() / 1e9 << endl;
    }
    int result = SUCCESS;

    for ( int num = 0; num < NITER; num++ )
    {
#pragma acc data present( P.P [0:2 * nxChunk * nyChunk * nzChunk * nChunk], this ) copy( result, err )
        //#pragma acc data  copy( result, err )
        {
            if ( MONITOR_MEM )
            {
                cout << "amount of free memory (0) = " << ( acc_get_free_memory() / 1e9 ) << endl;
            }
            if ( INITANALYTIC )
            {
#pragma acc parallel num_gangs( 1024 ) vector_length( VECLENGTH )
                initializeTrigonometric();
            }

#if ( POSS )
            if ( MONITOR_MEM )
            {
                cout << "amount of free memory after init = " << acc_get_free_memory() / 1e9 << endl;
            }
            if ( bc[0] != 'P' && bc[2] != 'P' )
            {
#pragma acc parallel num_gangs( 1024 ) vector_length( VECLENGTH )
                modifyRhsDirichlet();
            }

            if ( MONITOR_MEM )
            {
                cout << "amount of free memory 1 " << acc_get_free_memory() / 1e9 << endl;
            }

            t1 = MPI_Wtime();

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
            if ( GPUAW == 0 )
            {
                P.moveHostToDevice();
            }
            if ( MONITOR_MEM )
            {
                cout << "amount of free memory 2 " << acc_get_free_memory() / 1e9 << endl;
            }
            // perform FFT in x direction
            // remember, need to rearrange and restore each time

#if ( FFTX )
// step 2) change location of the array such that FFT can be performed on a contegeous array
#if ( COMM_PATTERN == 2 )
#pragma acc parallel num_gangs( trsps_gang0 ) vector_length( VECLENGTH )
            changeLocationXOverlap();
#else
#pragma acc parallel num_gangs( trsps_gang0 ) vector_length( VECLENGTH )
            changeLocationX();
#endif

            if ( MONITOR_MEM )
            {
                cout << "amount of free memory 3 " << acc_get_free_memory() / 1e9 << endl;
            }
// step 3) perform  FFT
#pragma acc parallel num_gangs( nSig0 ) vector_length( VECLENGTH )
            preprocessSignalAccordingly( 0, 0 );

            if ( MONITOR_MEM )
            {
                cout << "amount of free memory 4 " << acc_get_free_memory() / 1e9 << endl;
            }
            performTransformXdir();

            if ( MONITOR_MEM )
            {
                cout << "amount of free memory 5 " << acc_get_free_memory() / 1e9 << endl;
            }
#pragma acc parallel num_gangs( nSig0 ) vector_length( VECLENGTH )
            postprocessSignalAccordingly( 0, 0 );

            if ( MONITOR_MEM )
            {
                cout << "amount of free memory 6 " << acc_get_free_memory() / 1e9 << endl;
            }
// step 4) restore the array to original status before FFT
//#pragma acc parallel
#pragma acc parallel num_gangs( trsps_gang0 ) vector_length( VECLENGTH )
            restoreLocationX();

#endif

            if ( MONITOR_MEM )
            {
                cout << "amount of free memory 7 " << acc_get_free_memory() / 1e9 << endl;
            }
#if ( FFTY )
            // step 5) pencils with n(1,0,0) is converted to pencil with n(0,1,0)
            if ( GPUAW2 == 0 )
            {
                P.moveDeviceToHost();
            }
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

            if ( GPUAW2 == 0 )
            {
                P.moveHostToDevice();
            }
            // step 6) swaps X and Y coordinates, now X is Y and nx is ny

#pragma acc parallel num_gangs( trsps_gang0 ) vector_length( VECLENGTH )
            rearrangeX2Y();

            // step 7) change location of the array such that FFT can be performed on a contegeous array in the transverse direction

            //#pragma acc parallel
#if ( COMM_PATTERN == 2 )
#pragma acc parallel num_gangs( trsps_gang1 ) vector_length( VECLENGTH )
            // changeLocationYOverlap();
            changeLocationY();
#else
#pragma acc parallel num_gangs( trsps_gang1 ) vector_length( VECLENGTH )
            changeLocationY();
#endif

#pragma acc parallel num_gangs( nSig1 ) vector_length( VECLENGTH )
            preprocessSignalAccordingly( 1, 1 );

            // step 8) perform  FFT transform can be performed
            performTransformYdir();

            if ( MONITOR_MEM )
            {
                cout << "amount of free memory 8 " << acc_get_free_memory() / 1e9 << endl;
            }

            // acc_clear_freelists();
#pragma acc parallel num_gangs( nSig1 ) vector_length( VECLENGTH )
            postprocessSignalAccordingly( 1, 1 );
            // M.printX( myfile );

            // step 9) restore the array to original status before FFT
            //#pragma acc parallel

#pragma acc parallel num_gangs( trsps_gang1 ) vector_length( VECLENGTH )
            restoreLocationY();

#endif

            if ( MONITOR_MEM )
            {
                cout << "amount of free memory 9 " << acc_get_free_memory() / 1e9 << endl;
            }
#if ( SOLVE )
            // step 10) pencils with n(0,1,0) is converted to pencil with n(0,0,1)

            //   M.changeOwnershipPairwiseExchangeYZ();
            if ( GPUAW == 0 )
            {
                P.moveDeviceToHost();
            }

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
            if ( GPUAW == 0 )
            {
                P.moveHostToDevice();
            }
            // step 11) Customized multiBlock Thomas and periodic Thomas (with Sherman-Morrison modification)
            //#pragma acc parallel firstprivate( result ) reduction( + : result )
            //            result = solve();

            // generating two streams to handle the task parallel section of the code

            if ( SOLUTIONMETHOD == 0 )
            {
#pragma acc parallel num_gangs( gangTri ) vector_length( 1 )
                solveThmBatch( 0 );
                /*
                #pragma acc parallel
                                solveThmBatch( 1 );
                */
                /*
                #pragma acc parallel num_gangs( 1 ) async( 2 )
                                solveMG();

                #pragma acc parallel num_gangs( 1 ) async( 3 )
                                solveMGC();

                //#pragma acc wait(2,3)
                #pragma acc wait( 2 )
                #pragma acc wait( 3 )
                */

                if ( bc[0] == 'P' || bc[2] == 'P' )
                {
#pragma acc parallel num_gangs( gangTri ) vector_length( 1 )
                    solveThmBatch( 1 );
                }
            }

            else if ( SOLUTIONMETHOD == 1 )
            {
                /*
                //#pragma acc parallel async(123)
                #pragma acc parallel
                                solveThm( 0 );

                #pragma acc parallel
                                solveThm( 1 );
                */

                /*
                #pragma acc parallel num_gangs(nxChunk)
                            solveThmBatch( 0 );
                */

#pragma acc parallel num_gangs( gangTri ) vector_length( VECLENGTH )
                solvePCR( 0 );

                if ( bc[0] == 'P' || bc[2] == 'P' )
                {
#pragma acc parallel num_gangs( gangTri ) vector_length( VECLENGTH )
                    solvePCR( 1 );
                }
            }
            else if ( SOLUTIONMETHOD == 2 )
            {
#pragma acc parallel num_gangs( gangTri ) vector_length( VECLENGTH )
                solveCRP( 0 );

                if ( bc[0] == 'P' || bc[2] == 'P' )
                {
#pragma acc parallel num_gangs( gangTri ) vector_length( VECLENGTH )
                    solveCRP( 1 );
                }
            }

            else
            {
                // adding CUSPARSE

#if ( 0 )
                for ( int j = 0; j < nxChunk; j++ )
                {
                    for ( int i = 0; i < nyChunk; i++ )
                    {
                        //#pragma acc parallel
                        //                        eig = getEigenVal( i, j );

#pragma acc parallel async( 4 )
                        setDiag( i, j );

#pragma acc parallel async( 5 )
                        fillInArrayContig( i, j, 0 );

#pragma acc wait( 4 )
#pragma acc wait( 5 )

                        triDiagCusparse( dl, ds, du, tmpMGReal );

#pragma acc parallel
                        fillInArrayBack( i, j, 0 );
                    }
                }
/*
                   for ( int j = 0; j < nxChunk; j++ )
                  {
                        for ( int i = 0; i < nyChunk; i++ )
                        {
                #pragma acc parallel
                         eig = getEigenVal( i, j );
                #pragma acc parallel async(4)
                         setDiag(i,j);
                #pragma acc parallel async(5)
                         fillInArrayContig( i, j, 1 );
                #pragma acc wait(4)
                #pragma acc wait(5)
                         triDiagCusparse(dl,ds,du,tmpMGImag);
                #pragma acc parallel
                         fillInArrayBack( i, j, 1 );
                 }
                 }
*/
#endif
            }
            //#pragma acc wait(123)
            //            #pragma acc parallel
            //                        solveThm(1);
            /*
            #pragma acc update self (result)
                     // M.printX( myfile );
            */

            if ( MONITOR_MEM )
            {
                cout << "amount of free memory 10 " << acc_get_free_memory() / 1e9 << endl;
            }
            if ( result != SUCCESS )
            {
                cout << "Exit Code: " << THOMAS_FAIL << endl;
                cout << BLUE << PittPackGetErrorEnum( THOMAS_FAIL ) << RESET << endl;
                exit( 1 );
            }

            // step 12) pencils with n(0,0,1) is converted to pencil with n(0,1,0)
            // M.changeOwnershipPairwiseExchangeYZ();
            if ( GPUAW == 0 )
            {
                P.moveDeviceToHost();
            }
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

            if ( GPUAW == 0 )
            {
                P.moveHostToDevice();
            }
#endif

            if ( MONITOR_MEM )
            {
                cout << "amount of free memory 11 " << acc_get_free_memory() / 1e9 << endl;
            }
#if ( IFFTY )
// step 13) prepre for contigeuous FFT
#pragma acc parallel num_gangs( trsps_gang1 ) vector_length( VECLENGTH )
            changeLocationY();

            // step 14) perform IFFT

#pragma acc parallel num_gangs( nSig1 ) vector_length( VECLENGTH )
            preprocessSignalAccordinglyReverse( 1, 1 );

            performInverseTransformYdir();

#pragma acc parallel num_gangs( nSig1 ) vector_length( VECLENGTH )
            postprocessSignalAccordinglyReverse( 1, 1 );
            // step 15) restore the array to original status before IFFT

#pragma acc parallel num_gangs( trsps_gang1 ) vector_length( VECLENGTH )
            restoreLocationY();

// step 16) swaps X and Y coordinates, now X is Y and nx is ny
// need to redefine the function ??????????????????????????????????????????????????????????????????????????
// now again switching to x-direction arrangement so printX is fine
#pragma acc parallel num_gangs( trsps_gang0 ) vector_length( VECLENGTH )
            rearrangeX2YInverse();

            // step 17) pencils with n(0,1,0) is converted to pencil with n(1,0,0)

            if ( GPUAW2 == 0 )
            {
                P.moveDeviceToHost();
            }

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
            if ( GPUAW2 == 0 )
            {
                P.moveHostToDevice();
            }
#endif

            if ( MONITOR_MEM )
            {
                cout << "amount of free memory 12 " << acc_get_free_memory() / 1e9 << " GB " << endl;
            }
#if ( IFFTX )
            // step 18) change location of the array such that IFFT can be performed on a contegeous array

#if ( COMM_PATTERN == 2 )
#pragma acc parallel num_gangs( trsps_gang0 ) vector_length( VECLENGTH )
            changeLocationXOverlap();
            //  changeLocationX();
#else
#pragma acc parallel num_gangs( trsps_gang0 ) vector_length( VECLENGTH )
            changeLocationX();

#endif
            // step 19) perform  IFFT

#pragma acc parallel num_gangs( nSig0 ) vector_length( VECLENGTH )
            preprocessSignalAccordinglyReverse( 0, 0 );

            performInverseTransformXdir();

#pragma acc parallel num_gangs( nSig0 ) vector_length( VECLENGTH )
            postprocessSignalAccordinglyReverse( 0, 0 );

// step 20) restore the array to original status before FFT
#pragma acc parallel num_gangs( trsps_gang0 ) vector_length( VECLENGTH )
            restoreLocationX();

#pragma acc parallel num_gangs( nSig0 ) vector_length( VECLENGTH )
            rescale();

            if ( MONITOR_MEM )
            {
                cout << "amount of free memory 13 " << acc_get_free_memory() / 1e9 << " GB  " << endl;
            }
            // return back to the original set-up
            //
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
            //        setCoords( Xbox, 2 );
            //#pragma acc update device()

            t2 = MPI_Wtime();
            deT += ( t2 - t1 );

#endif

            if ( INCLUDE_ERROE_CAL_IN_TIMING == 1 )
            {
#pragma acc parallel vector_length( 32 ) reduction( max : err )
                err = getError();

                // cout<<" ****** "<<err<<endl;
            }

            if ( MONITOR_MEM )
            {
                cout << "amount of free memeory 14 " << acc_get_free_memory() / 1e9 << endl;
            }
#endif
        }
    }

    //#pragma acc update self(err)

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

    MPI_Reduce( &deT, &runTime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD );

    if ( myRank == 0 )
    {
        runTime = runTime / (double)NITER;
        runInfo();
        if ( PROFILE_COMM )
        {
            printf( "change ownership time =%lf percent of solution and take %lf seconds \n", t2_com / deT * 100., t2_com );
        }
    }

    if ( JIC )
    {
#if ( PITTPACKACC )
#pragma acc data present( P, R )
#endif
        {
#if ( PITTPACKACC )
#pragma acc parallel
#endif
            debug();
        }
    }

    P.moveDeviceToHost();

#if ( DEBUG2 )
    printX( myfile );
    myfile.close();
#endif
    cout << " *********************************" << endl;
#if ( DEBUG2 )
    //#if (1 )
    for ( int k = 0; k < 1; k++ )
    {
        for ( int j = 0; j < nyChunk; j++ )
        {
            for ( int i = 0; i < nxChunk; i++ )
            {
                cout << "\t " << P( i, j, k );
            }
            cout << endl;
        }
    }
#endif

#endif
}

#endif
