#include "pittpack.h"

void   test_cases( int NXCHUNK, int argcs, char *pArgs[] );
double test_case_DDDDDD( std::unique_ptr<PencilDcmp> &M );
double test_case_NNNNDD( std::unique_ptr<PencilDcmp> &M );
double test_case_NNNNNN( std::unique_ptr<PencilDcmp> &M );
double test_case_PPPPPP( std::unique_ptr<PencilDcmp> &M );
double test_case_NNPPPP( std::unique_ptr<PencilDcmp> &M,double *rhs );
double test_case_PPNNPP( std::unique_ptr<PencilDcmp> &M,double *rhs );
double test_case_DDPPPP( std::unique_ptr<PencilDcmp> &M,double *rhs );
double test_case_PPDDPP( std::unique_ptr<PencilDcmp> &M,double *rhs );
void   fillTrigonometricCosX(int myRank, int p0, double *X, int nx, int ny, int nz, double *rhs );
void   fillTrigonometricCosY(int myRank, int p0, double *X, int nx, int ny, int nz, double *rhs );
void   fillTrigonometricSinX(int myRank, int p0, double *X, int nx, int ny, int nz, double *rhs );
void   fillTrigonometricSinY(int myRank, int p0, double *X, int nx, int ny, int nz, double *rhs );

void test_cases( int NXCHUNK, int argcs, char *pArgs[] )
{
    int my_rank, com_size;
    int NYCHUNK = NXCHUNK;
    int NZCHUNK = NXCHUNK;

    MPI_Comm_rank( MPI_COMM_WORLD, &my_rank );
    MPI_Comm_size( MPI_COMM_WORLD, &com_size );

    int p0 = sqrt( com_size );
    int p1 = p0;

    if ( my_rank == 0 )
        cout << " p0 " << p0 << " " << p1 << endl;

    int Nx = NXCHUNK * p0;
    int Ny = NYCHUNK * p0;
    int Nz = NZCHUNK * p0;

    auto   M    = make_Poisson( argcs, pArgs, Nx, Ny, Nz );
    double X[6] = {0, 1.0, 0, 1.0, 0, 1.0};

    M->setBox( X );
    int dir = 2;
    M->setCoords( dir );

    int ncases = 20;

    double *err = new double[2 * ncases];

#if 0
    err[0] = test_case_NNNNDD( M );
    err[1] = test_case_DDDDDD( M );
    err[2] = test_case_NNNNNN( M );
    err[3] = test_case_PPPPPP( M );
    Nx = 2 * Nx;
    Ny = 2 * Ny;
    Nz = 2 * Nz;

    auto M0 = make_Poisson( argcs, pArgs, Nx, Ny, Nz );
    M0->setBox( X );
    M0->setCoords( dir );

    err[4] = test_case_NNNNDD( M0 );
    err[5] = test_case_DDDDDD( M0 );
    err[6] = test_case_NNNNNN( M0 );
    err[7] = test_case_PPPPPP( M0 );

    Nx = 2 * Nx;
    Ny = 2 * Ny;
    Nz = 2 * Nz;

    auto M1 = make_Poisson( argcs, pArgs, Nx, Ny, Nz );
    M1->setBox( X );
    M1->setCoords( dir );

    err[8]  = test_case_NNNNDD( M1 );
    err[9]  = test_case_DDDDDD( M1 );
    err[10] = test_case_NNNNNN( M1 );
    err[11] = test_case_PPPPPP( M1 );

    if ( my_rank == 0 )
    {
        std::cout << err[0] << std::endl;
        std::cout << err[1] << std::endl;
        std::cout << err[2] << std::endl;
        std::cout << err[3] << std::endl;
        std::cout << err[4] << std::endl;
        std::cout << err[5] << std::endl;
        std::cout << err[6] << std::endl;
        std::cout << err[7] << std::endl;
        std::cout << err[8] << std::endl;
        std::cout << err[9] << std::endl;
        std::cout << err[10] << std::endl;
        std::cout << err[11] << std::endl;
    }
    std::cout << "===============================================\n" << std::endl;
    std::cout << "----------- order of accuracy  ----------------\n" << std::endl;
    std::cout << "===============================================\n" << std::endl;
    std::cout << " NNNNDD "
              << " first refinement " << err[0] / err[4] << " second refinement " << err[4] / err[8] << std::endl;
    std::cout << " DDDDDD "
              << " first refinement " << err[1] / err[5] << " second refinement " << err[5] / err[9] << std::endl;
    std::cout << " NNNNNN "
              << " first refinement " << err[2] / err[6] << " second refinement " << err[6] / err[10] << std::endl;
    std::cout << " PPPPPP "
              << " first refinement " << err[3] / err[7] << " second refinement " << err[7] / err[11] << std::endl;

//==================================================================================
#endif
    cout<<" =============================== "<<endl;
    cout<<" ----------- 1D cases ---------- "<<endl;
    cout<<" =============================== "<<endl;
    Nx = NXCHUNK * p0;
    Ny = NYCHUNK * p0;
    Nz = NZCHUNK * p0;
    auto M2 = make_Poisson( argcs, pArgs, Nx, Ny, Nz );
    M2->setBox( X );
    M2->setCoords( dir );
    double *rhs=new double[Nx*Ny*Nz];

    fillTrigonometricCosX(my_rank, p0,X,Nx,Ny, Nz, rhs );
    err[0]= test_case_NNPPPP(M2,rhs);
    cout<< RED <<"err = "<<err[0]<<RESET<<endl;

    fillTrigonometricCosY(my_rank, p0,X,Nx,Ny, Nz, rhs );
    err[1]= test_case_PPNNPP(M2,rhs);
    cout<<RED<<" err = "<<err[1]<<RESET<<endl;

    fillTrigonometricSinX(my_rank,p0,X,Nx,Ny, Nz, rhs );
    err[2]= test_case_DDPPPP(M2,rhs);
    cout<<RED<<" err = "<<err[2]<<RESET<<endl;

    fillTrigonometricSinY(my_rank,p0,X,Ny,Nz, Nz, rhs );
    err[3]= test_case_PPDDPP(M2,rhs);
    cout<<RESET<<" err = "<<err[3]<<RESET<<endl;

    Nx = 2*Nx;
    Ny = 2*Ny;
    Nz = 2*Nz;

    auto M3 = make_Poisson( argcs, pArgs, Nx, Ny, Nz );
    M3->setBox( X );
    M3->setCoords( dir );
    
    double *rhs0=new double[Nx*Ny*Nz];
    fillTrigonometricCosX(my_rank, p0,X,Nx,Ny, Nz, rhs0 );
    err[4]= test_case_NNPPPP(M3,rhs0);
    cout<< RED <<"err = "<<err[4]<<RESET<<endl;

    fillTrigonometricCosY(my_rank, p0,X,Nx,Ny, Nz, rhs0 );
    err[5]= test_case_PPNNPP(M3,rhs0);
    cout<<RED<<" err = "<<err[5]<<RESET<<endl;

    fillTrigonometricSinX(my_rank,p0,X,Nx,Ny, Nz, rhs0 );
    err[6]= test_case_DDPPPP(M3,rhs0);
    cout<<RED<<" err = "<<err[6]<<RESET<<endl;

    fillTrigonometricSinY(my_rank,p0,X,Nx,Ny, Nz, rhs0 );
    err[7]= test_case_PPDDPP(M3,rhs0);
    cout<<RESET<<" err = "<<err[7]<<RESET<<endl;

if(my_rank==0)
{
    std::cout << "===============================================\n" << std::endl;
    std::cout << "----------- order of accuracy  ----------------\n" << std::endl;
    std::cout << "===============================================\n" << std::endl;
    std::cout << " NNPPPP "  << " first refinement " << err[0] / err[4] << std::endl;
    std::cout << " PPNNPP "  << " first refinement " << err[1] / err[5] <<  std::endl;
    std::cout << " DDPPPP "  << " first refinement " << err[2] / err[6] <<  std::endl;
    std::cout << " PPDDPP "  << " first refinement " << err[3] / err[7] <<  std::endl;
}

    delete[] err;
    delete[] rhs;
    delete[] rhs0;
}

double test_case_NNPPPP( std::unique_ptr<PencilDcmp> &M, double *rhs )
{
#if DEBUG
    std::cout << "=============================" << std::endl;
    std::cout << "           NNPPPP            " << std::endl;
    std::cout << "=============================" << std::endl;
#endif
    char mybc[6] = {'N', 'N', 'P', 'P', 'P', 'P'};
    M->assignBoundary( mybc );
    M->assignRhs(rhs);
    M->pittPack(); 
    double err=M->getErrorCos(0); 

return(err);
}

double test_case_PPNNPP( std::unique_ptr<PencilDcmp> &M, double *rhs )
{
#if DEBUG
    std::cout << "=============================" << std::endl;
    std::cout << "           PPNNPP            " << std::endl;
    std::cout << "=============================" << std::endl;
#endif
    char mybc[6] = {'P', 'P', 'N', 'N', 'P', 'P'};
    M->assignBoundary( mybc );
    M->assignRhs(rhs);
    M->pittPack(); 
    double err=M->getErrorCos(1); 
return(err);
}

double test_case_DDPPPP( std::unique_ptr<PencilDcmp> &M, double *rhs )
{
#if DEBUG
    std::cout << "=============================" << std::endl;
    std::cout << "           DDPPPP            " << std::endl;
    std::cout << "=============================" << std::endl;
#endif
    char mybc[6] = {'D', 'D', 'P', 'P', 'P', 'P'};
    M->assignBoundary( mybc );
    M->assignRhs(rhs);
    M->pittPack(); 
    double err=M->getErrorSin(0); 

return(err);
}

double test_case_PPDDPP( std::unique_ptr<PencilDcmp> &M, double *rhs )
{
#if DEBUG
    std::cout << "=============================" << std::endl;
    std::cout << "           PPNNPP            " << std::endl;
    std::cout << "=============================" << std::endl;
#endif
    char mybc[6] = {'P', 'P', 'D', 'D', 'P', 'P'};
    M->assignBoundary( mybc );
    M->assignRhs(rhs);
    M->pittPack(); 
    double err=M->getErrorSin(1); 
return(err);
}


double test_case_NNNNDD( std::unique_ptr<PencilDcmp> &M )
{
#if DEBUG
    std::cout << "=============================" << std::endl;
    std::cout << "           NNNNDD            " << std::endl;
    std::cout << "=============================" << std::endl;
#endif
    char mybc[6] = {'N', 'N', 'N', 'N', 'D', 'D'};
    M->assignBoundary( mybc );
    M->pittPack();
    return ( M->getError() );
}

double test_case_DDDDDD( std::unique_ptr<PencilDcmp> &M )
{
#if DEBUG
    std::cout << "=============================" << std::endl;
    std::cout << "           DDDDDD            " << std::endl;
    std::cout << "=============================" << std::endl;
#endif
    char mybc[6] = {'D', 'D', 'D', 'D', 'D', 'D'};
    M->assignBoundary( mybc );
    M->pittPack();
    return ( M->getError() );
}

double test_case_NNNNNN( std::unique_ptr<PencilDcmp> &M )
{
#if DEBUG
    std::cout << "=============================" << std::endl;
    std::cout << "           NNNNNN            " << std::endl;
    std::cout << "=============================" << std::endl;
#endif
    char mybc[6] = {'N', 'N', 'N', 'N', 'N', 'N'};
    M->assignBoundary( mybc );
    M->pittPack();
    return ( M->getError() );
}

double test_case_PPPPPP( std::unique_ptr<PencilDcmp> &M )
{
#if DEBUG
    std::cout << "=============================" << std::endl;
    std::cout << "           PPPPP            " << std::endl;
    std::cout << "=============================" << std::endl;
#endif
    char mybc[6] = {'P', 'P', 'P', 'P', 'P', 'P'};
    M->assignBoundary( mybc );
    M->pittPack();
    return ( M->getError() );
}

void fillTrigonometricCosX(int myRank, int p0, double *X, int Nx, int Ny, int nz, double *rhs )
{
    double pi = 4. * arctan( 1.0 );
    double x, y, z;
    double c1, c2, c3;

    int Nz = nz;

    double dxyz[3];
    dxyz[0] = -( X[0] - X[1] ) / ( Nx );
    dxyz[1] = -( X[2] - X[3] ) / ( Ny );
    dxyz[2] = -( X[4] - X[5] ) / ( Nz );

    double dx = ( X[1] - X[0] ) / double( p0 );
    double dy = ( X[3] - X[2] ) / double( p0 );

    double Xa = X[0] + ( myRank % p0 ) * dx;
    double Ya = X[2] + ( myRank / p0 ) * dy;
    double Za = X[4];


    c1 = dxyz[0];
    c2 = dxyz[1];
    c3 = dxyz[2];

    double shift = 1;

    for ( int k = 0; k < Nz; k++ )
    {
        z = Za + k * c3 + shift * c3 * 0.5;

        for ( int j = 0; j < Ny/p0; j++ )
        {
            y = Ya + j * c2 + shift * c2 * 0.5;

            for ( int i = 0; i < Nx/p0; i++ )
            {
                x = Xa + i * c1 + shift * c1 * .5;

                rhs[i + j * Nx/p0 + Nx/p0 * Ny/p0 * k] = -4. * pi * pi * cosine( 2. * pi * x );
            }
        }
    }
}

void fillTrigonometricCosY( int myRank, int p0, double *X, int Nx, int Ny, int nz, double *rhs )
{
    double pi = 4. * arctan( 1.0 );
    double x, y, z;
    double c1, c2, c3;

    //int Nx = nxChunk;
    //int Ny = nyChunk;
    int Nz = nz;

    double dxyz[3];
    dxyz[0] = -( X[0] - X[1] ) / ( Nx );
    dxyz[1] = -( X[2] - X[3] ) / ( Ny );
    dxyz[2] = -( X[4] - X[5] ) / ( Nz );
    c1      = dxyz[0];
    c2      = dxyz[1];
    c3      = dxyz[2];

    double dx = ( X[1] - X[0] ) / double( p0 );
    double dy = ( X[3] - X[2] ) / double( p0 );

    double Xa = X[0] + ( myRank % p0 ) * dx;
    double Ya = X[2] + ( myRank / p0 ) * dy;
    double Za = X[4];

    double shift = 1;

    for ( int k = 0; k < Nz; k++ )
    {
        z = Za + k * c3 + shift * c3 * 0.5;

        for ( int j = 0; j < Ny/p0; j++ )
        {
            y = Ya + j * c2 + shift * c2 * 0.5;

            for ( int i = 0; i < Nx/p0; i++ )
            {
                x = Xa + i * c1 + shift * c1 * .5;

                rhs[i + j * Nx/p0 + Nx/p0 * Ny/p0 * k] = -4. * pi * pi * cosine( 2. * pi * y );
            }
        }
    }
}

void fillTrigonometricSinX( int myRank, int p0, double *X, int Nx, int Ny, int nz, double *rhs )
{
    double pi = 4. * arctan( 1.0 );
    double x, y, z;
    double c1, c2, c3;

    int Nz = nz;

    double dxyz[3];
    dxyz[0] = -( X[0] - X[1] ) / ( Nx );
    dxyz[1] = -( X[2] - X[3] ) / ( Ny );
    dxyz[2] = -( X[4] - X[5] ) / ( Nz );

    double dx = ( X[1] - X[0] ) / double( p0 );
    double dy = ( X[3] - X[2] ) / double( p0 );

    double Xa = X[0] + ( myRank % p0 ) * dx;
    double Ya = X[2] + ( myRank / p0 ) * dy;
    double Za = X[4];

    c1 = dxyz[0];
    c2 = dxyz[1];
    c3 = dxyz[2];


    double shift = 1;

    for ( int k = 0; k < Nz; k++ )
    {
        z = Za + k * c3 + shift * c3 * 0.5;

        for ( int j = 0; j < Ny/p0; j++ )
        {
            y = Ya + j * c2 + shift * c2 * 0.5;

            for ( int i = 0; i < Nx/p0; i++ )
            {
                x = Xa + i * c1 + shift * c1 * .5;

                rhs[i + j * Nx/p0 + Nx/p0 * Ny/p0 * k] = -4. * pi * pi * sine( 2. * pi * x );
            }
        }
    }
}

void fillTrigonometricSinY( int myRank, int p0, double *X, int Nx, int Ny, int nz, double *rhs )
{
    double pi = 4. * arctan( 1.0 );
    double x, y, z;
    double c1, c2, c3;

    int Nz = nz;

    double dxyz[3];
    dxyz[0] = -( X[0] - X[1] ) / ( Nx );
    dxyz[1] = -( X[2] - X[3] ) / ( Ny );
    dxyz[2] = -( X[4] - X[5] ) / ( Nz );
    c1      = dxyz[0];
    c2      = dxyz[1];
    c3      = dxyz[2];

    double dx = ( X[1] - X[0] ) / double( p0 );
    double dy = ( X[3] - X[2] ) / double( p0 );

    double Xa = X[0] + ( myRank % p0 ) * dx;
    double Ya = X[2] + ( myRank / p0 ) * dy;
    double Za = X[4];

 

    double shift = 1;

    for ( int k = 0; k < Nz; k++ )
    {
        z = Za + k * c3 + shift * c3 * 0.5;

        for ( int j = 0; j < Ny/p0; j++ )
        {
            y = Ya + j * c2 + shift * c2 * 0.5;

            for ( int i = 0; i < Nx/p0; i++ )
            {
                x = Xa + i * c1 + shift * c1 * .5;

                rhs[i + j * Nx/p0 + Nx/p0 * Ny/p0 * k] = -4. * pi * pi * sine( 2. * pi * y );
            }
        }
    }
}
