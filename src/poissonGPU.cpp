#include "definitions.h"
#include "pencilDcmp.h"

#if ( PITTPACKACC )
void PoissonGPU::performTransformXdir() /*!< Is called on Host and Runs on GPU*/
{
    cufftHandle plan=NULL;
    double *ptr = P.P;

#pragma acc host_data use_device( ptr )
    {
        cufftPlan1d( &plan, nx, CUFFT_Z2Z, nyChunk * nzChunk );
       if(CUFFT_STREAMS){
        void *stream = acc_get_cuda_stream(acc_async_sync);
        cufftSetStream(plan, (cudaStream_t)stream);}
       if(CUFFT_SUCCESS!= cufftExecZ2Z( plan, (cufftDoubleComplex *)( ptr ), (cufftDoubleComplex *)( ptr ), CUFFT_FORWARD ))
       {
        cout << " Exit Code : " << PittPackGetErrorEnum( CUFFT_FAIL_X ) << endl; 
       }
        cufftDestroy( plan );
    }

}

void PoissonGPU::performInverseTransformXdir() /*!< Called on Host and Ran on GPU*/
{
    cufftHandle plan=NULL;
    double *ptr = P.P;
#pragma acc host_data use_device( ptr )
    {
        cufftPlan1d( &plan, nx, CUFFT_Z2Z, nyChunk * nzChunk );
       if(CUFFT_STREAMS){
        void *stream = acc_get_cuda_stream(acc_async_sync);
        cufftSetStream(plan, (cudaStream_t)stream);}
        if(CUFFT_SUCCESS!=cufftExecZ2Z( plan, (cufftDoubleComplex *)( ptr ), (cufftDoubleComplex *)( ptr ), CUFFT_INVERSE ))
        {
        cout << " Exit Code : " << PittPackGetErrorEnum( CUFFT_FAIL_INV_X ) << endl; 
        }
        cufftDestroy( plan );
    }
}
void PoissonGPU::performTransformYdir() /*!< Called on Host and Ran on GPU*/
{
    cufftHandle plan=NULL;
    double *ptr = P.P;
#pragma acc host_data use_device( ptr )
    {
        cufftPlan1d( &plan, ny, CUFFT_Z2Z, nxChunk * nzChunk );
       if(CUFFT_STREAMS)
       {
        void *stream = acc_get_cuda_stream(acc_async_sync);
        cufftSetStream(plan, (cudaStream_t)stream);
        }
       if(CUFFT_SUCCESS!= cufftExecZ2Z( plan, (cufftDoubleComplex *)( ptr ), (cufftDoubleComplex *)( ptr ), CUFFT_FORWARD ))
        {       
        cout << " Exit Code : " << PittPackGetErrorEnum( CUFFT_FAIL_Y ) << endl; 
        }
        cufftDestroy( plan );
    }
}

void PoissonGPU::performInverseTransformYdir() /*!< Called on Host and Ran on GPU*/
{
    cufftHandle plan=NULL;
    double *ptr = P.P;
#pragma acc host_data use_device( ptr )
    {
        cufftPlan1d( &plan, ny, CUFFT_Z2Z, nxChunk * nzChunk );
       if(CUFFT_STREAMS)
       {
        void *stream = acc_get_cuda_stream(acc_async_sync);
        cufftSetStream(plan, (cudaStream_t)stream);
        }
        if(CUFFT_SUCCESS!=cufftExecZ2Z( plan, (cufftDoubleComplex *)( ptr ), (cufftDoubleComplex *)( ptr ), CUFFT_INVERSE ))
        {
        cout << " Exit Code : " << PittPackGetErrorEnum( CUFFT_FAIL_INV_Y ) << endl; 
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
if(PIVOT==0)
{
        cusparseDgtsv_nopivot( handle, nz, 1, dl, ds, du, rhs, nz );
}
else
{
        cusparseDgtsv(handle, nz,1, dl, ds, du, rhs, nz);
}
    }

    cusparseDestroy( handle );
}

// const std::string currentDateTime() ;

void PoissonGPU::pittPack() /*!<called on CPU runs on GPU */
{
    //  struct timeval start_time, stop_time, elapsed_time;
    //  gettimeofday( &start_time, NULL );
    double t1;
    MPI_Barrier( MPI_COMM_WORLD );
    t1 = MPI_Wtime();

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
   
   long int initMem=acc_get_memory()/1e9; 
  
   cout<<"amount of available memory "<<initMem<<endl;

    double err = 0.0;
    finalErr = 0.0;
    double eig;
//    P.moveHostToDevice();
#if ( 1 )

    int trsps_gang0=MIN(50,iaxSize);
    int trsps_gang1=MIN(50, iaySize);

   int nSig0=MIN(50,nxChunk); 
   int nSig1=MIN(50,nyChunk); 
   initMem=acc_get_free_memory()/1e9;
    
    cout<<"amount of free memory before data region "<<acc_get_free_memory()/1e9<<endl;

    int result = SUCCESS;

    for ( int num = 0; num <1 ; num++ )
    {
#pragma acc data present( P.P[0 : 2 * nxChunk *nyChunk *nzChunk *nChunk], this) copy( result, err )
//#pragma acc data  copy( result, err )
        {

            cout<<"amount of free memory (0) = "<<(acc_get_free_memory()/1e9)<<endl;

#pragma acc parallel num_gangs(50)  vector_length(32)
            initializeTrigonometric();

#if ( POSS )
            cout<<"amount of free memory after init = "<<acc_get_free_memory()/1e9<<endl;
            if ( bc[0] != 'P' && bc[2] != 'P' )
            {
#pragma acc parallel num_gangs(50)  vector_length(VECLENGTH)
                modifyRhsDirichlet();
            }

            cout<<"amount of free memory 1 "<<acc_get_free_memory()/1e9<<endl;

            changeOwnershipPairwiseExchangeZX();

            if ( GPUAW == 0 )
            {
                P.moveHostToDevice();
            }
            cout<<"amount of free memory 2 "<<acc_get_free_memory()/1e9<<endl;
// perform FFT in x direction
// remember, need to rearrange and restore each time

#if ( FFTX )
// step 2) change location of the array such that FFT can be performed on a contegeous array
#pragma acc parallel num_gangs(trsps_gang0)  vector_length(VECLENGTH)
            changeLocationX();

            cout<<"amount of free memory 3 "<<acc_get_free_memory()/1e9<<endl;
// step 3) perform  FFT

#pragma acc parallel num_gangs(nSig0)  vector_length(VECLENGTH)
            preprocessSignalAccordingly( 0, 0 );

            cout<<"amount of free memory 4 "<<acc_get_free_memory()/1e9<<endl;

            performTransformXdir();

            cout<<"amount of free memory 5 "<<acc_get_free_memory()/1e9<<endl;

#pragma acc parallel num_gangs(nSig0)  vector_length(VECLENGTH)
            postprocessSignalAccordingly( 0, 0 );
            
            cout<<"amount of free memory 6 "<<acc_get_free_memory()/1e9<<endl;
// step 4) restore the array to original status before FFT
//#pragma acc parallel
#pragma acc parallel num_gangs(trsps_gang0)  vector_length(VECLENGTH)
            restoreLocationX();

#endif

            cout<<"amount of free memory 7 "<<acc_get_free_memory()/1e9<<endl;
#if ( FFTY )
            // step 5) pencils with n(1,0,0) is converted to pencil with n(0,1,0)
            if ( GPUAW2 == 0 )
            {
                P.moveDeviceToHost();
            }
            changeOwnershipPairwiseExchangeXY();
            if ( GPUAW2 == 0 )
            {
                P.moveHostToDevice();
            }
// step 6) swaps X and Y coordinates, now X is Y and nx is ny

#pragma acc parallel  num_gangs(trsps_gang0)  vector_length(VECLENGTH)
            rearrangeX2Y();

// step 7) change location of the array such that FFT can be performed on a contegeous array in the transverse direction

//#pragma acc parallel

#pragma acc parallel num_gangs(trsps_gang1)  vector_length(VECLENGTH)
            changeLocationY();

#pragma acc parallel  num_gangs(nSig1)  vector_length(VECLENGTH)
            preprocessSignalAccordingly( 1, 1 );

            // step 8) perform  FFT transform can be performed
            performTransformYdir();
            cout<<"amount of free memory 8 "<<acc_get_free_memory()/1e9<<endl;
           // acc_clear_freelists();
#pragma acc parallel  num_gangs(nSig1)  vector_length(VECLENGTH)
            postprocessSignalAccordingly( 1, 1 );
// M.printX( myfile );

// step 9) restore the array to original status before FFT
//#pragma acc parallel 

#pragma acc parallel num_gangs(trsps_gang1)  vector_length(VECLENGTH)
            restoreLocationY();

#endif

            cout<<"amount of free memory 9 "<<acc_get_free_memory()/1e9<<endl;
#if ( SOLVE )
            // step 10) pencils with n(0,1,0) is converted to pencil with n(0,0,1)

            //   M.changeOwnershipPairwiseExchangeYZ();
            if ( GPUAW == 0 )
            {
                P.moveDeviceToHost();
            }

            changeOwnershipPairwiseExchangeZX();

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

#pragma acc parallel num_gangs(gangTri) vector_length(1) 
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

#pragma acc parallel num_gangs(gangTri) vector_length(VECLENGTH)
            solvePCR( 0 );

            }
            else if(SOLUTIONMETHOD == 2 )
            {

#pragma acc parallel num_gangs(gangTri) vector_length(VECLENGTH) 
            solveCRP( 0 );

            }

            else
            {
                // adding CUSPARSE

#if(0)             
                for ( int j = 0; j < nxChunk; j++ )
                {
                    for ( int i = 0; i < nyChunk; i++ )
                    {

//#pragma acc parallel
//                        eig = getEigenVal( i, j );

#pragma acc parallel async( 4 )
                        setDiag( i,j );

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

            cout<<"amount of free memory 10 "<<acc_get_free_memory()/1e9<<endl;
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
            changeOwnershipPairwiseExchangeZX();

            if ( GPUAW == 0 )
            {
                P.moveHostToDevice();
            }
#endif

            cout<<"amount of free memory 11 "<<acc_get_free_memory()/1e9<<endl;
#if ( IFFTY )
// step 13) prepre for contigeuous FFT
#pragma acc parallel num_gangs(trsps_gang1)  vector_length(VECLENGTH)
            changeLocationY();

// step 14) perform IFFT

#pragma acc parallel  num_gangs(nSig1)  vector_length(VECLENGTH)
            preprocessSignalAccordinglyReverse( 1, 1 );

            performInverseTransformYdir();
#pragma acc parallel  num_gangs(nSig1)  vector_length(VECLENGTH)
            postprocessSignalAccordinglyReverse( 1, 1 );
// step 15) restore the array to original status before IFFT

#pragma acc parallel num_gangs(trsps_gang1)  vector_length(VECLENGTH)
            restoreLocationY();

// step 16) swaps X and Y coordinates, now X is Y and nx is ny
// need to redefine the function ??????????????????????????????????????????????????????????????????????????
// now again switching to x-direction arrangement so printX is fine
#pragma acc parallel num_gangs(trsps_gang0)  vector_length(VECLENGTH)
            rearrangeX2YInverse();

            // step 17) pencils with n(0,1,0) is converted to pencil with n(1,0,0)

            if ( GPUAW2 == 0 )
            {
                P.moveDeviceToHost();
            }
            changeOwnershipPairwiseExchangeXY();

            if ( GPUAW2 == 0 )
            {
                P.moveHostToDevice();
            }
#endif

            cout<<"amount of free memory 12 "<<acc_get_free_memory()/1e9<<" GB "<<endl;
#if ( IFFTX )
// step 18) change location of the array such that IFFT can be performed on a contegeous array

#pragma acc parallel num_gangs(trsps_gang0)  vector_length(VECLENGTH)
            changeLocationX();

// step 19) perform  IFFT

#pragma acc parallel  num_gangs(nSig0)  vector_length(VECLENGTH)
            preprocessSignalAccordinglyReverse( 0, 0 );

            performInverseTransformXdir();

#pragma acc parallel  num_gangs(nSig0)  vector_length(VECLENGTH)
            postprocessSignalAccordinglyReverse( 0, 0 );

// step 20) restore the array to original status before FFT
#pragma acc parallel num_gangs(trsps_gang0) vector_length(VECLENGTH)
            restoreLocationX();

#pragma acc parallel num_gangs(nSig0)  vector_length(VECLENGTH)
            rescale();

            cout<<"amount of free memory 13 "<<acc_get_free_memory()/1e9<<" GB  "<<endl;
            // return back to the original set-up
            changeOwnershipPairwiseExchangeZX();

//        setCoords( Xbox, 2 );
//#pragma acc update device()

#endif

            if ( INCLUDE_ERROE_CAL_IN_TIMING == 1 )
            {
#pragma acc parallel vector_length(50) reduction(max:err)
                err = getError();

                // cout<<" ****** "<<err<<endl;
            }

            cout<<"amount of free memeory 14 "<<acc_get_free_memory()/1e9<<endl;
#endif
        }
    }

    //#pragma acc update self(err)

    MPI_Allreduce( &err, &finalErr, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );

    // cout << YELLOW << "Error  (" << myRank << ") =" << err << " " << finalErr << RESET << endl;
    cout << YELLOW << "Error Per processor"
         << " " << finalErr / p0 / p0 << RESET << endl;

    finalErr = finalErr / p0 / p0;

    //  gettimeofday( &stop_time, NULL );
    //  timersub( &stop_time, &start_time, &elapsed_time ); // Unix time subtract routine

    double t2;
    MPI_Barrier( MPI_COMM_WORLD );
    t2 = MPI_Wtime();

    if ( myRank == 0 )
    {
        runTime = t2 - t1;
        runInfo();
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


