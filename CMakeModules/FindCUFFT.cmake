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
#   FFTW_USE_STATIC_LIBS    ... if true, only static libraries are found
#   FFTW_ROOT               ... if set, the libraries are exclusively searched
#                               under this path
#   FFTW_LIBRARY            ... fftw library to use
#   FFTW_INCLUDE_DIR        ... fftw include directory
#


#If environment variable FFTWDIR is specified, it has same effect as FFTW_ROOT

#if( NOT FFTW_ROOT AND $ENV{FFTWDIR} )
#  set( FFTW_ROOT $ENV{FFTW_DIR} )
#endif()

set( CUFFTW_ROOT $ENV{FFT_DIR} )
message("CUFFTW  :" ${CUFFTW_ROOT})


# Check if we can use PkgConfig
find_package(PkgConfig)

#Determine from PKG
if( PKG_CONFIG_FOUND AND NOT CUFFTW_ROOT )
  pkg_check_modules( PKG_FFTW QUIET "cufft" )
endif()

#Check whether to search static or dynamic libs
set( CMAKE_FIND_LIBRARY_SUFFIXES_SAV ${CMAKE_FIND_LIBRARY_SUFFIXES} )

message("Suffixes"  ${CMAKE_FIND_LIBRARY_SUFFIXES})

if( CUFFTW_ROOT )

 find_library(
    FFT_LIB
    NAMES "libcufft_static.a"
    PATHS ${CUFFTW_ROOT}
    PATH_SUFFIXES "lib" "lib64"
    NO_DEFAULT_PATH
)

  find_library(
    FFT_LIB
 #   NAMES "libcudart.so.7.5.18"
    PATHS ${CUFFTW_ROOT}/lib
    PATH_SUFFIXES "lib" "lib64"
    NO_DEFAULT_PATH
  )

  #find includes
  find_path(
    FFT_INCLUDES
    NAMES "cufftw.h"
    NAMES "cufft.h"
    PATHS ${CUFFTW_ROOT}/include
    PATH_SUFFIXES "include"
    NO_DEFAULT_PATH
  )

endif( CUFFTW_ROOT )

set(CUFFT_LIBRARIES ${FFT_LIB} )
set(CUFFT_INCLUDES ${FFT_INCLUDES} )

set( CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES_SAV} )

message("CUFFTW LIB  :" ${FFT_LIB})
message("CUFFTW INCLUDES  :" ${FFT_INCLUDES})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FFT DEFAULT_MSG CUFFT_INCLUDES CUFFT_LIBRARIES)

mark_as_advanced(CUFFT_INCLUDES CUFFT_LIBRARIES)

