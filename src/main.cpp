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

void test_cases(int NXCHUNK,  int argcs, char *pArgs[]) ;
double test_case_DDDDDD(std::unique_ptr<PencilDcmp> &M );
double test_case_NNNNDD(std::unique_ptr<PencilDcmp> &M );
double test_case_NNNNNN(std::unique_ptr<PencilDcmp> &M );
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

#if 0
    int my_rank, com_size;
    MPI_Comm_rank( MPI_COMM_WORLD, &my_rank );
    MPI_Comm_size( MPI_COMM_WORLD, &com_size );

    int p0 = sqrt( com_size );
    int p1 = p0;

    cout << " p0 " << p0 << " " << p1 << endl;

    int Nx = NXCHUNK * p0;
    int Ny = NYCHUNK * p0;
    int Nz = NZCHUNK * p0;

     auto M=make_Poisson( argcs, pArgs, Nx, Ny, Nz );
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
    M->assignBoundary( mybc );

    double a[3] = {0, 0, 0};
//    cout << RED << " myRank " << my_rank << " a[3]= " << a[0] << " " << a[1] << RESET << endl;

    // double X[6] = {-1.0, 1.0, -1.0, 1.0, -1.0, 1.0};
    double X[6] = {0, 1.0, 0, 1.0, 0, 1.0};

    M->setBox( X );

   int dir = 2;

    cout << "xxxxxxxxxxxxx set coords xxxxxxxxxxxxxxxxxxxxx" << endl;
    M->setCoords( dir );

   // step 1) pencils with n(0,0,1) is converted to pencil with n(1,0,0)

    double *rhs = nullptr;

#if ( 1 )
#if ( !EXACT )
//#if ( OPENACC )

        if ( INITANALYTIC == 0 )
        {
            double *rhs = new double[Nx * Ny * Nz];
            M->fillTrigonometric( rhs );
                M->assignRhs( rhs );
        }
        M->pittPack();

#endif
// final write should be in the dir=2 as we transform all the data back to its original situation
// if you dont call poisson, if you do set it equal to dir=0

#if ( 1 )
    cout << " for writing out" << endl;
    M->setCoords( dir );
    if ( I_O == 1 )
    {
        M->IO( 1, dir, 0 );
    }

   if ( INCLUDE_ERROE_CAL_IN_TIMING == 0 )
    {
        cout << RED << " err = " << M->getError() << RESET << endl;
    }
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
#endif 

//   test_case_NNNNDD(NXCHUNK,argcs,pArgs );
//   test_case_DDDDDD(NXCHUNK,argcs,pArgs );


    test_cases(NXCHUNK,  argcs, pArgs );

    return ( 0 );
};


void test_cases(int NXCHUNK,  int argcs, char *pArgs[] )
{
    int my_rank,com_size;
    int NYCHUNK = NXCHUNK;
    int NZCHUNK = NXCHUNK;
  
    MPI_Comm_rank( MPI_COMM_WORLD, &my_rank );
    MPI_Comm_size( MPI_COMM_WORLD, &com_size );

    int p0 = sqrt( com_size );
    int p1 = p0;

    cout << " p0 " << p0 << " " << p1 << endl;

    int Nx = NXCHUNK * p0;
    int Ny = NYCHUNK * p0;
    int Nz = NZCHUNK * p0;

    auto M=make_Poisson( argcs, pArgs, Nx, Ny, Nz );
    double a[3] = {0, 0, 0};
    double X[6] = {0, 1.0, 0, 1.0, 0, 1.0};
 
    M->setBox( X );
    int dir = 2;
    M->setCoords( dir );

   int ncases=9;

   double *err=new double[2*ncases];
   err[0]=test_case_NNNNDD(M);
   err[1]=test_case_DDDDDD(M);
   err[2]=test_case_NNNNNN(M);

    Nx = 2*Nx;
    Ny = 2*Ny;
    Nz = 2*Nz;

    auto M0=make_Poisson( argcs, pArgs, Nx, Ny, Nz );
    M0->setBox( X );
    M0->setCoords( dir );
 
   err[3]=test_case_NNNNDD(M0);
   err[4]=test_case_DDDDDD(M0);
   err[5]=test_case_NNNNNN(M0);

    Nx = 2*Nx;
    Ny = 2*Ny;
    Nz = 2*Nz;

    auto M1=make_Poisson( argcs, pArgs, Nx, Ny, Nz );
    M1->setBox( X );
    M1->setCoords( dir );
 
   err[6]=test_case_NNNNDD(M1);
   err[7]=test_case_DDDDDD(M1);
   err[8]=test_case_NNNNNN(M1);


   std::cout<<err[0]<<std::endl;
   std::cout<<err[1]<<std::endl;
   std::cout<<err[2]<<std::endl;
   std::cout<<err[3]<<std::endl;
   std::cout<<err[4]<<std::endl;
   std::cout<<err[5]<<std::endl;
   std::cout<<err[6]<<std::endl;
   std::cout<<err[7]<<std::endl;
   std::cout<<err[8]<<std::endl;
//   std::cout<<test_case_DDDDDD(M)<<std::endl;

}

double test_case_NNNNDD(std::unique_ptr<PencilDcmp> &M )
{   
    std::cout<<"============================="<<std::endl;
    std::cout<<"           NNNNDD            "<<std::endl;
    std::cout<<"============================="<<std::endl;
    char mybc[6] = {'N', 'N', 'N', 'N', 'D', 'D'};
    M->assignBoundary( mybc );
    M->pittPack();
    return(M->getError()); 
}


double test_case_DDDDDD(std::unique_ptr<PencilDcmp> &M )
{
    std::cout<<"============================="<<std::endl;
    std::cout<<"           DDDDDD            "<<std::endl;
    std::cout<<"============================="<<std::endl;
    char mybc[6] = {'D', 'D', 'D', 'D', 'D', 'D'};
    M->assignBoundary( mybc );
    M->pittPack();
    return(M->getError()); 
}


double test_case_NNNNNN(std::unique_ptr<PencilDcmp> &M )
{
    std::cout<<"============================="<<std::endl;
    std::cout<<"           NNNNNN            "<<std::endl;
    std::cout<<"============================="<<std::endl;
    char mybc[6] = {'N', 'N', 'N', 'N', 'N', 'N'};
    M->assignBoundary( mybc );
    M->pittPack();
    return(M->getError()); 
}












