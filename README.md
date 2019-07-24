# PittPack
PittPack is an Open-Source Poisson’s Equation Solver for Extreme-Scale Computing with Accelerators
It is developed at University of Pittsburgh as part of the GEM3D project funded by NSF

# Automatically configures the code to run on Stampede2 (TACC), Comet (SD) and Bridges (PSC)  


## On the Linux Terminal do:
### source config.sh 
## When prompted respond by entering 0 or 1    
###  0 will configure PittPack for CPU clusters 
###  1 will configure PittPack for GPU clusters


for complete documentation visit www.pittpack.com


## List of Required Libraries
### cmake 
### FFTW3 
### PGI 
### HDF5
### MPI 

## To Build 
### run ./config.sh
### cd build
### cmake ..
### make -j 4
### executable will be placed in /bin

## To run
### mpirun -np numberofProcs ./bin/PittPack nx ny nz 
 
