# PittPack
PittPack is an Open-Source Poisson’s Equation Solver for Extreme-Scale Computing with Accelerators <br/>
It is developed as part of the GEM3D project funded by NSF at University of Pittsburgh, PA, USA. 

# Configuration 
automatically configures the code for Stampede2 (TACC), Comet (SD) and Bridges (PSC) clusters 

## On the Linux Terminal do:
source config.sh <br/> 
When prompted respond by entering 0 or 1    
* 0 : will configure PittPack for CPU clusters 
* 1 : will configure PittPack for GPU clusters



# Installation
PittPack requires the following libraries
  * cmake 
  * FFTW3 
  * PGI 
  * HDF5
  * MPI 

#  Build  
Perform the following steps
  * source config.sh
  * cd build
  * cmake ..
  * make -j 4
  * executable will be placed in /bin

## To run
mpirun -np numberofProcs ./bin/PittPack nx ny nz 
 
# Detailed usage
for complete documentation visit www.pittpack.com
