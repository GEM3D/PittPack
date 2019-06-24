#ifndef _DEFINITIONS_H_
#define _DEFINITIONS_H_
#include "fftw3.h"
#include "hdf5.h"
#include "params.h"
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <mpi.h>
#include <sstream>
#include <stdexcept>
#include <string>
#include <sys/time.h>
#include <time.h>
#include <vector>

#if ( PITTPACKACC )
#include "cusparse.h"
#include "math.h"
#include "openacc.h"
#include <cublas_v2.h>
#include <cuda_runtime.h>
#endif

typedef int integer;

typedef unsigned int unit;
typedef double       PittPackReal;

using namespace std;

#if ( SHORT_ )
typedef unsigned short sint;
#else
// typedef unsigned int sint;
typedef unsigned long sint;
#endif

template <class T>
std::string FormatWithCommas( T value )
{
    std::stringstream ss;
    ss.imbue( std::locale( "" ) );
    ss << std::fixed << value;
    return ss.str();
}

#define RED "\033[01;31m"
#define GREEN "\033[22;32m"
#define YELLOW "\033[22;33m"
#define BLUE "\033[22;34m"
#define MAGENTA "\033[22;35m"
#define CYAN "\033[22;36m"
#define RESET "\033[22;0m"

#endif
