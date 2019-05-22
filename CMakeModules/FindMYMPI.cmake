# taken from https://github.com/libigl/eigen/blob/master/cmake/FindMPI.cmake
#
# - Find the MPI library
# Usage:
#   find_package(MPI [REQUIRED] [QUIET] )
#     
# It sets the following variables:
#   MPI_FOUND               ... true if fftw is found on the system
#   MPI_LIBRARIES           ... full path to fftw library
#   MPI_INCLUDES            ... fftw include directory
#
# The following variables will be checked by the function
#   MPI_USE_STATIC_LIBS    ... if true, only static libraries are found
#   MPI_ROOT               ... if set, the libraries are exclusively searched
#                               under this path
#   MPI_LIBRARY            ... fftw library to use
#   MPI_INCLUDE_DIR        ... fftw include directory
#

set( MPI_ROOT $ENV{MPI_ROOT} )

if(NOT MPI_ROOT)
	set(MPI_ROOT $ENV{MPIHOME})
endif()

if(NOT MPI_ROOT)
	set(MPI_ROOT $ENV{MPI_DIR})
endif()

if(NOT MPI_ROOT)
	set(MPI_ROOT $ENV{MPI_HOME})
endif()

# Stampede (TACC)
if(NOT MPI_ROOT)
	set(MPI_ROOT $ENV{I_MPI_ROOT})
endif()

# Titan (ORNL)
if(NOT MPI_ROOT)
	set(MPI_ROOT $ENV{MPICH_DIR})
endif()


message("MPI" ${MPI_ROOT})

# Check if we can use PkgConfig
find_package(PkgConfig)

#Determine from PKG
if( PKG_CONFIG_FOUND AND NOT MPI_ROOT )
	pkg_check_modules( PKG_MPI QUIET "mpi" )
endif()

#Check whether to search static or dynamic libs
set( CMAKE_FIND_LIBRARY_SUFFIXES_SAV ${CMAKE_FIND_LIBRARY_SUFFIXES} )

message("Suffixes"  ${CMAKE_FIND_LIBRARY_SUFFIXES})

find_library(
	MPI_LIB
	NAMES "libmpi.so" "libmpich.so"
	PATHS ${MPI_ROOT}
	PATH_SUFFIXES "lib" "lib64"
	NO_DEFAULT_PATH
)
#find includes
find_path(
	MPI_INCLUDES
	NAMES "mpi.h"
	PATHS ${MPI_ROOT}/include
	PATH_SUFFIXES "include"
	NO_DEFAULT_PATH
)


set(MPI_LIBRARIES ${MPI_LIB})

if(MPIL_LIB)
	set(MPI_LIBRARIES ${MPI_LIBRARIES} ${MPIL_LIB})
endif()

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(MPI DEFAULT_MSG MPI_INCLUDES MPI_LIBRARIES)

mark_as_advanced(MPI_INCLUDES MPI_LIBRARIES )

