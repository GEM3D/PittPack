# Synopsis
**PittPack** is an Open-Source Poisson’s Equation Solver for Extreme-Scale Computing with Accelerators <br/>
The main goal is to solve the Poisson's equation with second order accuracy on directionally uniform Cartesian gird for conventional
as well as accelerated clusters. It uses two FFT transforms in x and y directions along with the tridiagonal solve in z direction.
Due to the FFT it limits the boundary conditions in the x and y directions to one of the following combinations
  * Periodic
  * Neumann-Neumann
  * Neumann-Dirichlet
  * Dirichlet-Dirichlet
 

## Features
  * Hybrid MPI/OpenACC parallelization
  * Chunked-Pencil Decomposition
  * Low-memory communication pattern option via pairwise-exchange algorithm
  * User-friendly interface   
 

## Configuration 
PittPack provides config.sh which can be used to automatically configure the code for Stampede2 (https://www.tacc.utexas.edu/systems/stampede2), Comet (https://www.sdsc.edu/services/hpc/hpc_systems.html) and Bridges (https://www.psc.edu/resources/computing/bridges) clusters 

## Linux 
```
source config.sh 
```
When prompted respond by entering 0 or 1    
* (0) : will configure PittPack for CPU clusters 
* (1) : will configure PittPack for GPU clusters



## Installation
PittPack requires the following libraries
  * cmake 
  * MPI 
  * FFTW3 
  * PGI  
  * HDF5 with Parallel IO
  * ParaView (or VisIt) for visualization of the result

##  Build  
PittPack uses *CMakeLists.txt* and *CMakeModules* folder to detect the library paths. <br/>
These two components are crucial for complilation of PittPack.
Perform the following steps
```
  cd build
  cmake ..
  make 
```
The executable will be placed in the /bin folder


## Run
```
mpirun -np N ./bin/PittPack nx ny nz 
```
  * N: Number of processes (squared integer)
  * nx: Number of elements in X-direction
  * ny: Number of elements in Y-direction
  * nz: Number of elements in Z-direction
 
## Visualization
  * The output is written to the /soln folder 
  * Paraview should be used to visualize the solution
  * Simply open the file ending with xmf in soln/ 


## Directory structure
```
PittPack
│   README.md
│   CMakeLists.txt    
│   LICENSE
│   config.sh
│
└─── CMakeModules
│   │   FindPGI.cmake: cmake script to find PGI 
│   │   FindMYMPI.cmake: cmake script to find MPI
│   │   FindMYHDF5.cmake: cmake script to find HDF5
│   │   FindGOOGLETEST.cmake: cmake script to find Google Test
│   │   FindFFTW.cmake: cmake script to find FFTW
│
└─── src
│   │   chunkedArray.cpp:  Abstracts access patterns 
│   │   communicate.cpp:   MPI communications
│   │   mathFunction.cpp:  Defines math functions for kernel generation by PGI
│   │   signalProc.cpp:    Performs FFT transforms
│   │   poissonCPU.cpp:    Inherits from class PencilDcmp and specialized for CPU
│   │   poissonGPU.cpp:    Inherits from class PencilDcmp and specialized for GPU
│   │   pencilDcmp.cpp:    Incorporates Decomposition strategy and communication patterns
│   │   triDiag.cpp:       Class for tridiagonal solvers      
│   │   phdf5.cpp:         Class for handling IO with hdf5     
│   │  
│   └─── include (headers)
│       │   chunkedArray.hpp
│       │   communicate.hpp  
│       │   mathFunction.hpp
│       │   signalProc.hpp
│       │   poissonCPU.hpp
│       │   poissonGPU.hpp
│       │   pencilDcmp.hpp
│       │   triDiag.hpp 
│       │   phdf5.hpp
│       │   params.h 
│       │   definitions.h 
│   
└─── bin
│       PittPack (executable)  
│  
│
└─── build   
│       Will be populated by cmake   
│  
│
└─── soln 
│   │   Pxdmf3d1.h5: outputs the file in hdf5 format 
│   │   Pxdmf3d1.xmf: meta data to be used by ParaView
│   
└─── archives
 
```

## Details
For complete documentation visit www.pittpack.com

## Notes 
We welcome any feedbacks by the users and developers <br/>
Please read the LICENSE file for how to use this software

## Acknowledgements
**PittPack** is developed as part of the GEM3D project funded by NSF at University of Pittsburgh, PA, USA. 


## Contributors
  * Jaber J. Hasbestan (jaber@pitt.edu)
  * Inanc Senocak (senocak@pitt.edu)

