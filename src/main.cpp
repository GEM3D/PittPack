#include "pittpack.h"

/*!  \mainpage
 *
 *  \details  "GPU accelerated FFT-based Poisson Solver" <br>
 * This project is part of the <STRONG>NSF-GEM3D </STRONG>Award No. 1440638. <br>
 *
 *
 *
 *             Required Libraries:
 *               CMAKE
 *               CUFFT
 *               FFTW3
 *               HDF5
 *               PGI (v.18.1)
 *               MPI
 *
 *
 *             Usage:
 *             mpirun -np NP  nx ny nz
 *
 * @authors     Jaber J. Hasbestan and Inanc Senocak (PI) <br>
 * @date        Jan 2019
 * \details
 * Copyright (c)  by the Authors 2019 \n
 * Swanson School of Engineering \n
 * Department of Mechanical Engineering and Materials Science \n
 * University of Pittsburgh \n
 *
 * Distribution is subjected to GLP 3.0 license
 *  \image html sample_poisson.png
 */

void checkCudaSupport();
void parse();

int main( int argcs, char *pArgs[] )
{
    if ( argcs < 3 )
    {
        cout << " Exit Code : " << PittPackGetErrorEnum( INPUT_ARGS ) << endl;
        exit( 0 );
    }

    for ( int i = 0; i < argcs; i++ )
    {
        cout << pArgs[i] << "\t";
    }
    cout << endl;

    int NXCHUNK = atoi( pArgs[1] );
    int NYCHUNK = atoi( pArgs[2] );
    int NZCHUNK = atoi( pArgs[3] );

    if ( MPI_Init( &argcs, &pArgs ) != MPI_SUCCESS )
    {
        cout << " Exit Code : " << PittPackGetErrorEnum( MPI_INIT_FAIL ) << endl;
        exit( 1 );
    }

    int my_rank, com_size;
    MPI_Comm_rank( MPI_COMM_WORLD, &my_rank );
    MPI_Comm_size( MPI_COMM_WORLD, &com_size );
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

#if ( OPENACC )
    //    #if ( 0 )
    if ( my_rank == 0 )
    {
        cout << " ====================================== " << endl;
        cout << " ====================================== " << endl;
        cout << " OPENACC is defined "
                "============"
             << endl;
        cout << " ====================================== " << endl;
        cout << " ====================================== " << endl;
    }
    PittPackResult result;
    result = OPENACC_Init( my_rank, com_size );

    if ( SUCCESS != result )
    {
        cout << RED << " Rank(" << my_rank << ") > Exit Code : " << result << RESET << endl;
        cout << BLUE << PittPackGetErrorEnum( result ) << RESET << endl;
        exit( 1 );
    }

#endif
    // HostToDeviceAssign(my_rank,com_size);

    int p0 = sqrt( com_size );
    int p1 = p0;

    cout << " p0 " << p0 << " " << p1 << endl;

    int Nx = NXCHUNK * p0;
    int Ny = NYCHUNK * p0;
    int Nz = NZCHUNK * p0;

#if ( OPENACC )
    cout << " GPU solving " << endl;
     auto M=make_Poisson( argcs, pArgs, Nx, Ny, Nz );

#else
    cout << " CPU solving "
         << " Nx " << Nx << " " << Ny << " " << Nz << " p0 " << p0 << endl;
    cout << RED << " Make Sure -ta=pinned is turned-off " << RESET << endl;
     auto M=make_Poisson( argcs, pArgs, Nx, Ny, Nz );
#endif

     // char mybc[6] = {'P', 'P', 'P', 'P', 'D', 'D'};
    // char mybc[6] = {'D', 'D', 'P', 'P', 'P', 'P'};
    //  char mybc[6] = {'P', 'P', 'P', 'P', 'D', 'D'};
    // char mybc[6] = {'P', 'P', 'P', 'P', 'P', 'P'};
    // it is illposed
    //  char mybc[6] = {'N', 'N', 'N', 'N', 'N', 'N'};
    // char mybc[6] = {'D', 'D', 'D', 'D', 'D', 'D'};

     char mybc[6] = {'N', 'N', 'N', 'N', 'D', 'D'};
    // char mybc[6] = {'D', 'D', 'D', 'D', 'P', 'P'};
    //char mybc[6] = {'P', 'P', 'P', 'P', 'N', 'N'};
    //       char mybc[6] = {'D', 'D', 'D', 'D', 'P', 'P'};
    //        char mybc[6] = {'D', 'D', 'D', 'D', 'N', 'N'};
    // char mybc[6] = {'N', 'N', 'N', 'N', 'N', 'N'};
    std::cout << mybc[0] << " " << mybc[1] << " " << mybc[2] << " " << mybc[3] << " " << mybc[4] << " " << mybc[5] << std::endl;
    //M.assignBoundary( mybc );
    M->assignBoundary( mybc );
    // testMpiClass(MPI_COMM_WORLD);

    double a[3] = {0, 0, 0};
    //      M.modifyEigForDirichlet(0,0,a);
    cout << RED << " myRank " << my_rank << " a[3]= " << a[0] << " " << a[1] << RESET << endl;

    cout << "============================ " << endl;

    // This is only for neighbrhood collectives
    M->constructConnectivity();

    M->graphCreate();

    // double X[6] = {-1.0, 1.0, -1.0, 1.0, -1.0, 1.0};
    double X[6] = {0, 1.0, 0, 1.0, 0, 1.0};

    M->setBox( X );

    // check IO for dir 0

#if ( DEBUG1 )
    ofstream myfile;

    std::string filename = "data";
    filename.append( to_string( my_rank ) );
    //  ofstream myfile;
    myfile.open( filename );
#endif
    // if(my_rank==0)
    // debugging for single block first
    int dir = 2;

    cout << "xxxxxxxxxxxxx set coords xxxxxxxxxxxxxxxxxxxxx" << endl;
    M->setCoords( dir );

#if ( DEBUG1 )
    myfile << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx" << endl;
    myfile << "     Start" << endl;

    //    if ( my_rank == this_rank )
    {
        M.printX( myfile );
    }
#endif
    // step 1) pencils with n(0,0,1) is converted to pencil with n(1,0,0)

    double *rhs = nullptr;

#if ( 1 )
#if ( !EXACT )
//#if ( OPENACC )

    //    for ( int i = 0; i < 10; i++ )
    {
        //        M.initializeTrigonometric();
        if ( INITANALYTIC == 0 )
        {
            double *rhs = new double[Nx * Ny * Nz];
            //M.fillTrigonometric( rhs );
            M->fillTrigonometric( rhs );
            /*
            for(int i=0;i<Nx*Ny*Nz;i++)
            {
            cout<<rhs[i]<<endl;
            }
             */
            //M.assignRhs( rhs );
            M->assignRhs( rhs );
        }
        // M.print();
        //M.pittPack();
        M->pittPack();
        //  delete [] rhs;
    }

#endif
// final write should be in the dir=2 as we transform all the data back to its original situation
// if you dont call poisson, if you do set it equal to dir=0
#if ( EXACT )
    dir = 2;
#else
    // used to be 1
    dir = 2;
#endif

#if ( 1 )
    cout << " for writing out" << endl;
    //M.setCoords( dir );
    M->setCoords( dir );
    if ( I_O == 1 )
    {
        //  M.changeOwnershipPairwiseExchangeZX();
        //  M.changeOwnershipPairwiseExchangeXY();
        //M.IO( 1, dir, 0 );
        M->IO( 1, dir, 0 );
    }

//   M.eigenVal( myfile );
#if ( DEBUG1 )
    myfile.close();
#endif
    if ( INCLUDE_ERROE_CAL_IN_TIMING == 0 )
    {
        //cout << RED << " err = " << M.getError() << RESET << endl;
        cout << RED << " err = " << M->getError() << RESET << endl;

    }
//    M.runInfo();
#endif

    #endif
    if ( INCLUDE_ERROE_CAL_IN_TIMING == 1 )
    {
        cout << RED << "Warning : output file show the error not the solution  " << RESET << endl;
    }
    if ( rhs != nullptr )
    {
        delete[] rhs;
    }

    return ( 0 );
};

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



