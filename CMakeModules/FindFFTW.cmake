# taken from https://github.com/libigl/eigen/blob/master/cmake/FindFFTW.cmake
# The part for finding suffixes is removed since it found the wrong sufixes
# - Find the FFTW library
#
# Usage:
#   find_package(FFTW [REQUIRED] [QUIET] )
#     
# It sets the following variables:
#   FFTW_FOUND               ... true if fftw is found on the system
#   FFTW_LIBRARIES           ... full path to fftw library
#   FFTW_INCLUDES            ... fftw include directory
#
# The following variables will be checked by the function
#   FFTW_LIBRARY            ... fftw library to use
#   FFTW_INCLUDE_DIR        ... fftw include directory
#

#  Different Clusters use different environment variables
#  here we include the most common ones used in clusters all over the country

if(NOT FFTW_ROOT)
	set( FFTW_ROOT $ENV{FFTWHOME} )
endif()

if(NOT FFTW_ROOT)
	set( FFTW_ROOT $ENV{FFTWROOT} )
endif()

if(NOT FFTW_ROOT)
	set( FFTW_ROOT $ENV{TACC_FFTW3_DIR} )
endif()

# Titan (ORNL)
if(NOT FFTW_ROOT)
	set( FFTW_ROOT "/opt/cray/fftw/3.3.4.11/interlagos/" )
endif()

if(NOT FFTW_ROOT)
	set( FFTW_ROOT $ENV{FFTW_DIR} )
endif()


message("FFTW" ${FFTW_ROOT})

# Check if we can use PkgConfig
find_package(PkgConfig)

#Determine from PKG
if( PKG_CONFIG_FOUND AND NOT FFTW_ROOT )
	pkg_check_modules( PKG_FFTW QUIET "fftw3" )
endif()

#Check whether to search static or dynamic libs
#set( CMAKE_FIND_LIBRARY_SUFFIXES_SAV ${CMAKE_FIND_LIBRARY_SUFFIXES} )

message("Suffixes"  ${CMAKE_FIND_LIBRARY_SUFFIXES})

if( FFTW_ROOT )

find_library(
	FFTW_LIB
	NAMES "libfftw3.so"
	PATHS ${FFTW_ROOT}
	PATH_SUFFIXES "lib" "lib64"
	NO_DEFAULT_PATH
)

#find includes
find_path(
	FFTW_INCLUDES
	NAMES "fftw3.h"
	PATHS ${FFTW_ROOT}/include
	PATH_SUFFIXES "include"
	NO_DEFAULT_PATH
)

else()

find_library(
	FFTW_LIB
	NAMES "libfftw3.a"
	PATHS ${PKG_FFTW_LIBRARY_DIRS} ${LIB_INSTALL_DIR}
)

find_path(
	FFTW_INCLUDES
	NAMES "fftw3.h"
	PATHS ${PKG_FFTW_INCLUDE_DIRS} ${INCLUDE_INSTALL_DIR}
)

endif( FFTW_ROOT )

set(FFTW_LIBRARIES ${FFTW_LIB} ${FFTWF_LIB})

if(FFTWL_LIB)
	set(FFTW_LIBRARIES ${FFTW_LIBRARIES} ${FFTWL_LIB})
endif()

set( CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES_SAV} )

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(FFTW DEFAULT_MSG FFTW_INCLUDES FFTW_LIBRARIES)

mark_as_advanced(FFTW_INCLUDES FFTW_LIBRARIES FFTW_LIB FFTWF_LIB FFTWL_LIB)

