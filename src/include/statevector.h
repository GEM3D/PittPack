#ifndef _STATEVECTOR_H_
#define _STATEVECTOR_H_
//#include "cufft.h"
#if ( 0 )
#include <algorithm>
#include <bitset>
#include <cmath>
#include <complex>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <functional>
#include <stack>
#include <stdexcept>
#include <unordered_map>
#include <vector>
#endif
#include <assert.h>

#include <iostream>

using namespace std;

class Fourier;

/*!
 *
 * \class Complex Pressure Overloads Some Operators and makes
 *  the padded complex use easy
 * This class initiates a cell centered pressure field
 *
 * */

class ComplexPressure
{
    private:
    int Nx;
    int Ny;
    int Nz;
    int arraySize;
    double *__restrict__ P;
    //    double *P;

    public:
    ComplexPressure(){};
    ComplexPressure( int n ); /*! Constructor  */
    void moveHostToDevice();  /*! Constructor  */
    void moveDeviceToHost();  /*! Constructor  */
    void allocate( int n );
    void allocate( int *n );
    int  size();
    void print();
    ~ComplexPressure();
#pragma acc routine
    double &operator()( int i, int j, int k, int index );
#pragma acc routine
    double &operator()( int i, int j, int k );

    friend class Fourier;
    friend void poisson();

    //    friend cufftResult cufftPlan1d(cufftHandle *plan, int rank, cufftType type, int batch);

    //    friend  cufftResult cufftPlanMany(cufftHandle *plan, int rank, int *n, int *inembed, int istride, int idist, int *onembed, int
    // ostride, int odist, cufftType type, int batch);

    void getAddress( double *rt );

    void performTransform();

    //#pragma acc routine
    //    int operator()( int i, int j, int k );
};

/*!
 *
 * \class Double  Pressure Overloads Some Operators and makes
 *  the padded complex use easy
 * This class initiates a cell centered pressure field
 *
 * */

class Pressure
{
    private:
    int Nx;
    int Ny;
    int Nz;
    int arraySize;
    double *__restrict__ P;

    public:
    Pressure(){};
    Pressure( int n );       /*! Constructor  */
    void moveHostToDevice(); /*! Constructor  */
    void moveDeviceToHost(); /*! Constructor  */
    void allocate( int n );
    void allocate( int *n );
    int  size();
    void print();
    ~Pressure();
#pragma acc routine
    double &operator()( int i, int j, int k );

    //#pragma acc routine
    //    int operator()( int i, int j, int k );
};

#endif
