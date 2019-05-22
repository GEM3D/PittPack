
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
	PATHS ${CUFFTW_ROOT}/lib
	PATH_SUFFIXES "lib" "lib64"
	NO_DEFAULT_PATH
  )

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

