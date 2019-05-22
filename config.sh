#!/bin/bash

GREEN='\033[22;32m'
RESET='\033[0m'
RED='\033[01;31m'

#module purge


string=$HOSTNAME
if [[ $string == *"comet"* ]]; then
    export cluster=comet
    echo $cluster
elif [[ $string == *"bridges"* ]]; then
    export cluster=bridges
    echo $cluster
elif [[ $string == *"stampede"* ]]; then
    export cluster=stampede
    echo $cluster
elif [[ $string == *"crc"* ]]; then
    export cluster=crc
    echo $cluster
elif [[ $string == *"titan"* ]]; then
    export cluster=titan
    echo $cluster
fi


printf "${RED}\n \n       System is =  $cluster ${RESET}\n \n"

read -p "To Configure for GPU press 1 else enter 0:  " GPU
echo $GPU

if [ $cluster == "bridges" ]; then

    #module purge
    export HDF5_ROOT=/pylon5/ac561fp/hasbesta/packages/hdf5-1.8.17/
    module load psc_path/1.1
    module load slurm/default
    module load xdusage/2.0-4
    module load cmake/3.7.2
    module load fftw3/3.3.4


    if [ $GPU -lt 1 ]; then
        module load /opt/packages/openmpi/openmpi-2.0.1/modulefiles/gcc_openmpi-2.0.1
        export CXX=/opt/packages/openmpi/openmpi-2.0.1/bin/mpiCC;
        export CC=/opt/packages/openmpi/openmpi-2.0.1/bin/mpicc;
        export GPU1=0;
        export _OPENACC=0;
    else
        module load mpi/pgi_openmpi.18.1
        export MPI_ROOT=/opt/packages/pgi/linux86-64/2018/mpi/openmpi-2.1.2/
        export CC=$MPI_ROOT/bin/mpicc;
        export CXX=$MPI_ROOT/bin/mpiCC;
        export GPU1=1;
        export _OPENACC=1;
    fi

fi


if [ $cluster == "comet" ]; then
    module purge
    export MODULEPATH=/opt/modulefiles/applications/.gnu:/opt/modulefiles/mpi/.gnu:/opt/modulefiles/mpi:/opt/modulefiles/compilers:/opt/modulefiles/applications:/usr/share/Modules/modulefiles:/etc/modulefiles:/share/apps/compute/modulefiles/mpi:

    export MODULEPATH=$MODULEPATH:/share/apps/compute/modulefiles:/share/apps/compute/modulefiles/mpi:/share/apps/compute/modulefiles/applications

    export MODULEPATH=$MODULEPATH:/share/apps/compute/doxygen/bin

    #export MODULEPATH=/opt/modulefiles/applications/.gnu:/opt/modulefiles/mpi/.gnu:/opt/modulefiles/mpi:/opt/modulefiles/compilers:/opt/modulefiles/applications:/usr/share/Modules/modulefiles:/etc/modulefiles:/share/apps/compute/modulefiles/mpi:

    #system needs
    module load cmake/3.9.1
    module load gnutools
    module load gnu/6.2.0
    #module load fftw/3.3.4

    #export HDF5HOME=/opt/hdf5/intel/openmpi_ib
    ##module load .intel/mvapich2_ib/2.1

    if [ $GPU -lt 1 ]; then
        module load openmpi_ib/2.1.1_gnu
        module load hdf5/1.10.2gnu 
        module load fftw
	#module load fftw/3.3.4 
	#export FFTWHOME=/opt/fftw/3.3.4/gnu/openmpi_ib;
        export CC=$MPIHOME/bin/mpicc
        export CXX=$MPIHOME/bin/mpiCC
        export GPU1=0;
        export _OPENACC=0;
    else


        export MODULEPATH=/share/apps/gpu/pgi/v18.4/modulefiles:$MODULEPATH
        module purge
        module load gnubase
        module load pgi/18.4
        module load openmpi/2.1.2
        module load hdf5/1.10.2
        module load fftw
	### Set FFTW manually
	#export FFTWHOME=/opt/fftw/3.3.4/pgi/openmpi_ib 
	#export PATH=/opt/fftw/3.3.4/pgi/openmpi_ib/bin:$PATH
	#export LD_LIBRARY_PATH=/opt/fftw/3.3.4/pgi/openmpi_ib/lib:$LD_LIBRARY_PATH
	#export MODULEPATH=/share/apps/compute/modulefiles:$MODULEPATH
        module load gnu/4.9.2
        module load cmake
        module cuda
	#module load fftw/3.3.4
	#cp /share/apps/gpu/pgi/v18.4/linux86-64/18.4/bin/localrc.comet.gcc620 ~/.mypgcpprc
        export CC=$MPIHOME/bin/mpicc;
        export CXX=$MPIHOME/bin/mpiCC;
        export GPU1=1;
        export _OPENACC=1;
    fi

fi



if [ $cluster == "stampede" ]; then

    module load intel/18.0.2  impi/18.0.2
    module load phdf5/1.8.16


    #export HDF5HOME=/opt/hdf5/intel/openmpi_ib
    ##module load .intel/mvapich2_ib/2.1

    if [ $GPU -lt 1 ]; then
	#module load openmpi_ib/2.0.4_gnu
        module load fftw3/3.3.8
        export CC=$TACC_IMPI_DIR/intel64/bin/mpicc
        export CXX=$TACC_IMPI_DIR/intel64/bin/mpicxx
        export GPU1=0;
    else
        echo "No GPU is available on Stampede"
    fi

fi


if [ $cluster == "crc" ]; then

     module purge
     module load intel/2018.2.199
     module load intel-mpi/2018.2.199
     module load gcc/6.3.0
     module load hdf5/1.10.0
     module load fftw/3.3.8
     module load  cmake/3.7.1



    if [ $GPU -lt 1 ]; then
        export CC=$I_MPI_ROOT/intel64/bin/mpicc
        export CXX=$I_MPI_ROOT/intel64/bin/mpiicpc
        export GPU1=0;
    else
        export GPU1=1
        module load gcc/5.4.0
        module load openmpi/2.0.2
        module load pgi/18.10
        export CC=/ihome/crc/install/pgi/18.10/linux86-64/18.10/mpi/openmpi-2.1.2/bin/mpicc
        export CXX=/ihome/crc/install/pgi/18.10/linux86-64/18.10/mpi/openmpi-2.1.2/bin/mpiCC
    fi

fi


if [ $cluster == "titan" ]; then

    #system needs
    module load cmake3
    module unload PrgEnv-pgi/5.2.82
    module unload PrgEnv-gnu/5.2.82

    export CRAYPE_LINK_TYPE=dynamic

    if [ $GPU -lt 1 ]; then
        module load PrgEnv-gnu/5.2.82
        module load cray-mpich/7.6.3
        module load cray-hdf5-parallel
        module load fftw
        export GPU1=0;
    else
        module load PrgEnv-pgi
        module swap pgi pgi/18.10.0
        module load cray-mpich/7.6.3
        module load cray-hdf5-parallel
        module load cudatoolkit
        export GPU1=1;
    fi

fi

printf "${RED} Creating Required Directories ${RESET}\n \n"

if [ ! -d "./soln" ]; then
    mkdir soln
fi

if [ ! -d "./build" ]; then
    mkdir build
fi

printf "${GREEN} Required modules have been loaded for your system! ${RESET}\n \n"
printf "${GREEN} Change directory to buid/ and type cmake .. and make sure it finds all packages ${RESET}\n \n"
printf "${GREEN} After successful cmake, type make to compile the code. make -j n to compile in parallel with n processes ${RESET}\n \n"


################################################################

#                  FORMATING the SHELL SCRIPT

################################################################

#1 Load the file into Emacs "emacs config.sh"
#2 Press Ctrl-space at the top of the file
#3 Move the cursor to the bottom of the file
#4 Press Alt-X and type untabify then return
#4 Press Alt-X and type indent-region then return


