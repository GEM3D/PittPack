#include "pittpack.h"


void test_cases(int NXCHUNK,  int argcs, char *pArgs[]) ;
double test_case_DDDDDD(std::unique_ptr<PencilDcmp> &M );
double test_case_NNNNDD(std::unique_ptr<PencilDcmp> &M );
double test_case_NNNNNN(std::unique_ptr<PencilDcmp> &M );
double test_case_PPPPPP(std::unique_ptr<PencilDcmp> &M );


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

   int ncases=20;

   double *err=new double[2*ncases];
   err[0]=test_case_NNNNDD(M);
   err[1]=test_case_DDDDDD(M);
   err[2]=test_case_NNNNNN(M);
   err[3]=test_case_PPPPPP(M);

    Nx = 2*Nx;
    Ny = 2*Ny;
    Nz = 2*Nz;

    auto M0=make_Poisson( argcs, pArgs, Nx, Ny, Nz );
    M0->setBox( X );
    M0->setCoords( dir );
 
   err[4]=test_case_NNNNDD(M0);
   err[5]=test_case_DDDDDD(M0);
   err[6]=test_case_NNNNNN(M0);
   err[7]=test_case_PPPPPP(M0);

    Nx = 2*Nx;
    Ny = 2*Ny;
    Nz = 2*Nz;

    auto M1=make_Poisson( argcs, pArgs, Nx, Ny, Nz );
    M1->setBox( X );
    M1->setCoords( dir );
 
   err[8]=test_case_NNNNDD(M1);
   err[9]=test_case_DDDDDD(M1);
   err[10]=test_case_NNNNNN(M1);
   err[11]=test_case_PPPPPP(M1);


   std::cout<<err[0]<<std::endl;
   std::cout<<err[1]<<std::endl;
   std::cout<<err[2]<<std::endl;
   std::cout<<err[3]<<std::endl;
   std::cout<<err[4]<<std::endl;
   std::cout<<err[5]<<std::endl;
   std::cout<<err[6]<<std::endl;
   std::cout<<err[7]<<std::endl;
   std::cout<<err[8]<<std::endl;
   std::cout<<err[9]<<std::endl;
   std::cout<<err[10]<<std::endl;
   std::cout<<err[11]<<std::endl;
   delete[] err;
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

double test_case_PPPPPP(std::unique_ptr<PencilDcmp> &M )
{
    std::cout<<"============================="<<std::endl;
    std::cout<<"           PPPPP            "<<std::endl;
    std::cout<<"============================="<<std::endl;
    char mybc[6] = {'P', 'P', 'P', 'P', 'P', 'P'};
    M->assignBoundary( mybc );
    M->pittPack();
    return(M->getError()); 
}














