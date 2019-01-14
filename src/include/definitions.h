#ifndef _DEFINITIONS_H_
#define _DEFINITIONS_H_
#include "fftw3.h"
#include "hdf5.h"
#include "params.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <iostream>
#include <list>
#include <memory>
#include <mpi.h>
#include <sstream>
#include <stack>
#include <stdexcept>
#include <string>
#include <sys/time.h>
#include <time.h>
#include <unistd.h>
#include <unordered_map>
#include <vector>
#ifdef _OPENACC
#include "math.h"
#include "openacc.h"
#endif

typedef int integer;

typedef unsigned int unit;
// typedef float real;
typedef double real;

using namespace std;

#if ( SHORT )
typedef unsigned short sint;
#else
// typedef unsigned int sint;
typedef int sint;
#endif


template<class T>
std::string FormatWithCommas(T value)
{
    std::stringstream ss;
    ss.imbue(std::locale(""));
    ss << std::fixed << value;
    return ss.str();
}



/*
using real = double;
using uint = unsigned int;
using integer = int;
*/

#define RED "\033[01;31m"
#define GREEN "\033[22;32m"
#define YELLOW "\033[22;33m"
#define BLUE "\033[22;34m"
#define MAGENTA "\033[22;35m"
#define CYAN "\033[22;36m"
#define RESET "\033[22;0m"

#endif
