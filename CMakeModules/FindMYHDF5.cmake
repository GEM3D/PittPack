

set(HDF5_ROOT $ENV{HDF5_ROOT})

if(NOT HDF5_ROOT)
	set(HDF5_ROOT $ENV{HDF5HOME})
endif()

if(NOT HDF5_ROOT)
	set(HDF5_ROOT $ENV{HDF5_DIR})
endif()

if(NOT HDF5_ROOT)
	set(HDF5_ROOT $ENV{TACC_HDF5_DIR})
endif()

message("HDF5_ROOT  :" ${HDF5_ROOT})

# Check if we can use PkgConfig
find_package(PkgConfig)

#Determine from PKG
if( PKG_CONFIG_FOUND AND NOT HDF5_ROOT )
	pkg_check_modules( PKG_FFTW QUIET "hdf5" )
endif()

#Check whether to search static or dynamic libs
set( CMAKE_FIND_LIBRARY_SUFFIXES_SAV ${CMAKE_FIND_LIBRARY_SUFFIXES} )

message("Suffixes"  ${CMAKE_FIND_LIBRARY_SUFFIXES})

if( HDF5_ROOT )

find_library(
	HDF5_LIB
	NAMES "libhdf5.so" "libhdf5-shared.so"
	PATHS ${HDF5_ROOT}
	PATH_SUFFIXES "lib" "lib64"
	NO_DEFAULT_PATH
)

#  find_library(
#    HDF5_LIB
#    PATHS ${HDF5_ROOT}/lib
#    NAMES "libhdf5.so"
#    PATH_SUFFIXES "lib" "lib64"
#    NO_DEFAULT_PATH
#  )

#find includes
find_path(
	HDF5_INCLUDES
	NAMES "hdf5.h"
	PATHS ${HDF5_ROOT}/include
	PATH_SUFFIXES "include"
	NO_DEFAULT_PATH
)

endif( HDF5_ROOT )


set(HDF5_hdf5_LIBRARY_RELEASE ${HDF5_LIB} )
set(HDF5_C_INCLUDE_DIR ${HDF5_INCLUDES} )
set( CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES_SAV} )
set(HDF5_LIBRARIES ${HDF5_LIB} )
set(HDF5_INCLUDE_DIRS ${HDF5_INCLUDES})

message("HDF5 LIB  :" ${HDF5_LIBRARIES})
message("HDF5 INCLUDES  :" ${HDF5_INCLUDES})

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(HDF5 DEFAULT_MSG HDF5_INCLUDES HDF5_LIBRARIES)

mark_as_advanced(HDF5_C_INCLUDE_DIR HDF5_hdf5_LIBRARY_REALEASE)



