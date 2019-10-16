#include "pittpack.h"
    
void sanityCheck(  int argcs, char *pArgs[]  )
{

// call test cases

    int NXCHUNK = 32;
    int NYCHUNK = NXCHUNK;
    int NZCHUNK = NXCHUNK;

    if ( MPI_Init( &argcs, &pArgs ) != MPI_SUCCESS )
    {
        cout << " Exit Code : " << PittPackGetErrorEnum( MPI_INIT_FAIL ) << endl;
        exit( 1 );
    }

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
     char mybc[6] = {'N', 'N', 'N', 'N', 'D', 'D'};
     std::cout << mybc[0] << " " << mybc[1] << " " << mybc[2] << " " << mybc[3] << " " << mybc[4] << " " << mybc[5] << std::endl;
     M->assignBoundary( mybc );
    double a[3] = {0, 0, 0};
    M->constructConnectivity();
    M->graphCreate();
    double X[6] = {0, 1.0, 0, 1.0, 0, 1.0};
    M->setBox( X );
    int dir = 2;
    M->setCoords( dir );
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
#endif


// final write should be in the dir=2 as we transform all the data back to its original situation
// if you dont call poisson, if you do set it equal to dir=0

    if ( rhs != nullptr )
    {
        delete[] rhs;
    }

}

void test_case_NNNNDD(PencilDcmp  *M,int NXCHUNK,int NYCHUNK,int NZCHUNK, int p0)
{

    int p1 = p0;


    int Nx = NXCHUNK * p0;
    int Ny = NYCHUNK * p0;
    int Nz = NZCHUNK * p0;

    char mybc[6] = {'N', 'N', 'N', 'N', 'D', 'D'};
     std::cout << mybc[0] << " " << mybc[1] << " " << mybc[2] << " " << mybc[3] << " " << mybc[4] << " " << mybc[5] << std::endl;

    M->assignBoundary( mybc );
    double a[3] = {0, 0, 0};
    double X[6] = {0, 1.0, 0, 1.0, 0, 1.0};
 //  M->setBox( X );
 //   int dir = 2;
 //   M->setCoords( dir );
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
#endif


// final write should be in the dir=2 as we transform all the data back to its original situation
// if you dont call poisson, if you do set it equal to dir=0

    if ( rhs != nullptr )
    {
        delete[] rhs;
    }

}



