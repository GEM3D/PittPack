cmake_minimum_required(VERSION 2.8.11)

project(pittPack)

#disable in source build
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)

message(STATUS "cmake is building on a ${CMAKE_SYSTEM_NAME} system")
set(dir ${PROJECT_SOURCE_DIR})
message("project source dir" ${dir})

set(CMAKE_FILES_DIRECTORY ${dir}/build)
message("Files Directory" ${CMAKE_FILES_DIRECTORY})

#define some colors (Unix-based systems)

if(NOT WIN32)
  string(ASCII 27 Esc)
  set(ColourReset "${Esc}[m")
  set(ColourBold  "${Esc}[1m")
  set(Red         "${Esc}[31m")
  set(Green       "${Esc}[32m")
  set(Yellow      "${Esc}[33m")
  set(Blue        "${Esc}[34m")
  set(Magenta     "${Esc}[35m")
  set(Cyan        "${Esc}[36m")
  set(White       "${Esc}[37m")
  set(BoldRed     "${Esc}[1;31m")
  set(BoldGreen   "${Esc}[1;32m")
  set(BoldYellow  "${Esc}[1;33m")
  set(BoldBlue    "${Esc}[1;34m")
  set(BoldMagenta "${Esc}[1;35m")
  set(BoldCyan    "${Esc}[1;36m")
  set(BoldWhite   "${Esc}[1;37m")
endif()


# get the number of physical procesors 

message("project binary dir" ${EXECUTABLE_OUTPUT_PATH})
message("CMAKE SOURCE DIR:" ${CMAKE_SOURCE_DIR})
message("CMAKE_BINARY_DIR:" ${CMAKE_BINARY_DIR})

# Out-of-Source build is required by amrGem, this if statement prevents
# in-source-builds 

if ( ${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
    message( FATAL_ERROR "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there. You may need to remove CMakeCache.txt." )
endif()


message("cmake Home Dir:" ${CMAKE_HOME_DIRECTORY})
message("cmake Build Dir:" ${CMAKE_BUILD_DIRECTORY})
message("cmake Binary Dir:" ${CMAKE_BINARY_DIR})
message("cmake source dir:" ${CMAKE_SOURCE_DIR})


############################################################################
#
#                    SECTION: FINDING-PACKAGE
#
#        1) HDF5 have native FindXXX.cmake's provided by cmake
#        2) For Zoltan and ParMETIS one needs to specify FindXXX.cmake 
#        3) These non-native to cmake modules are included in CMakeModules
#
############################################################################

# specify the path for Modules 
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/CMakeModules/")

#find HDF5
#find_package(HDF5 ${HDF5_MIN_VERSION} REQUIRED)
#if(HDF5_FOUND)
#include_directories(${HDF5_INCLUDE_DIRS})
#message( "${Red}HDF5 Library Found: ${ColourReset}" ${HDF5_INCLUDE_DIRS})
#message("${Green} HDF5 Library Found: ${ColourReset}")
#endif()

find_package(MYMPI ${MPI_MIN_VERSION} REQUIRED)
if(MPI_FOUND)
include_directories(${MPI_INCLUDE_PATH})
message( "${Red} MPI Include Path : ${ColourReset}" ${MPI_INCLUDES})
message( "${Red} MPI Library Path : ${ColourReset}"  ${MPI_LIBRARIES})
#message("${Red} MPI Library Found : ${ColourReset}")
endif()


#find FFTW
find_package(FFTW REQUIRED)
if(FFTW_FOUND)
include_directories(${FFTW_INCLUDES})
message("${Red} FFTW Library Found: ${ColourReset}")
endif()
message("${Red} ${FFTW_INCLUDES} ${ColourReset}")

#find HDF5
find_package(MYHDF5 ${HDF5_MIN_VERSION} REQUIRED)
if(HDF5_FOUND)
include_directories(${HDF5_INCLUDE_DIRS})
message( "${Red} HDF5 Library Found: ${ColourReset}" ${HDF5_INCLUDE_DIRS})
message("${Green} 2. HDF5 Library : Found ${ColourReset}")
endif()


#find_package(PGI REQUIRED)
#if(PGI_FOUND)
#include_directories(${PGI_INCLUDES})
#message("${Red} PGI Library Found: ${ColourReset}")
#endif()
#message("${Red} ${PGI_INCLUDES} ${ColourReset}")

#find_package(CUFFT REQUIRED)
#if(CUFFT_FOUND)
#include_directories(${CUFFT_INCLUDES})
#message("${Red}  CUFFT Found: ${ColourReset}")
#endif()

#message("${Green} CUFFT Library Found: ${CUFFT_LIBRARIES} ${ColourReset}")
#message("${Green} CUFFT includes:  ${CUFFT_INCLUDES}  ${ColourReset}")

###########################################################################
#
#                    SECTION: BUILD AND INSTALL
#
#        1) Build all the *.cpp and *.c files in "/src" directory
#        2) Link with external package libraries 
#        3) Install the final executable in "/bin"
#
############################################################################

#OPTION(GPU "Option description" $ENV{GPU1})
set(GPU $ENV{GPU1})
if(GPU)
add_definitions(-DOPENACC=1)
endif(GPU)
message("${Red} GPU is ${GPU} ${ColourReset}" )

# This is performed by source modules.sh
#if(GPU)
#set(CC /opt/packages/pgi/linux86-64/2018/mpi/openmpi-2.1.2/bin/mpicc)
#set(CXX /opt/packages/pgi/linux86-64/2018/mpi/openmpi-2.1.2/bin/mpiCC)
#else(GPU)
#set(CXX /opt/packages/openmpi/openmpi-2.0.1/bin/mpiCC)
#set(CC /opt/packages/openmpi/openmpi-2.0.1/bin/mpicc)
#endif(GPU)


message("${Green} CXX ${CXX}   ${ColourReset} ")


include_directories(${dir}/src/include)
# dont forget to modify the address of the source files
#set(CMAKE_BUILD_FILES_DIRECTORY ${dir})
#set(CMAKE_BUILD_DIRECTORY ${dir}/bin)      

file(GLOB OBJ "${dir}/src/*.cpp" )
set(EXECUTABLE_OUTPUT_PATH ${dir}/bin)

# Define all header files for linking, this is required such that one doesnt
# have to hardcode the directories for I/O operations and hence ensure
# portability

#set(HEADER_FILES "${dir}/src/include/tridiag.h" "${dir}/src/include/statevector.h" "${dir}/src/include/fft.h" )
set(HEADER_FILES    "${dir}/src/include/communicate.h"  "${dir}/src/include/pencilDcmp.h"  "${dir}/src/include/definitions.h"  "${dir}/src/include/datatype.h" "${dir}/src/include/chunkedArray.h" "${dir}/src/include/phdf5.h"   "${dir}/src/include/params.h"  "${dir}/src/include/signalProc.h"  "${dir}/src/include/mathFunction.h"  "${dir}/src/include/multiGrid.h" )

#set(HEADER_FILES  "${dir}/src/include/*")



message("header files" ${HEADER_FILES})

# insert -pg instead of -g for checking with gprof
# gprof bin/amrGem > out.text  
# -Minfo=severe turns off the warnings
message("GPU" $ENV{GPU1})


if($ENV{GPU1})
## for GPU

 SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DOMPI_SKIP_MPICXX -acc -fast -Minfo=accel -Mlarge_arrays -Mcudalib=cufft,cusparse  -ta=tesla,gvmode   -std=c++11")
 SET(CMAKE_EXEC_LINKER_FLAGS "${CMAKE_EXEC_LINKER_FLAGS} -DOMPI_SKIP_MPICXX -acc -fast -Minfo=accel -Mlarge_arrays  -ta=tesla,gvmode -Mcudalib=cufft,cusparse -std=c++11" )

else()
   SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++1y -O3  -Wno-unknown-pragmas  ")
   #SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++1y -O3  -Wall ")
# SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -acc -O3 -Minfo=accel -Mcudalib=cufft  -ta=tesla:cuda8.0,pinned   -Minfo=ccff -Minform=severe -std=c++11")
#  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -acc -O3 -Minfo=accel -Mcudalib=cufft  -ta=tesla:cuda9.0   -Minfo=ccff -Minform=severe -std=c++11")
#  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -acc   -Minfo=accel -Mcudalib=cufft  -ta=tesla:cuda8.0   -Minfo=ccff  -std=c++11")
#    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -acc -g  -Minfo=accel -Mcudalib=cufft  -ta=multicore   -Minfo=ccff -Minform=severe -std=c++11")
#    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -Mcudalib=cufft  std=c++0x")   -ta=tesla:managed -
#   SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -std=c++1y -pg -O3 -Wno-unused-variable -Wno-unknown-pragmas -Wno-sign-compare ")
#    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=c++1y  -O3 -Wunused-variable ")
#    SET(CMAKE_EXEC_LINKER_FLAGS "${CMAKE_EXEC_LINKER_FLAGS} -Wall -std=gnu++0x" )
    SET(CMAKE_EXEC_LINKER_FLAGS "${CMAKE_EXEC_LINKER_FLAGS} -O3 -std=c++1y -Wno-unknown-pragmas -Wno-sign-compare " )
#   SET(CMAKE_EXEC_LINKER_FLAGS "${CMAKE_EXEC_LINKER_FLAGS}   -Minfo=accel -ta=tesla,cuda9.0,pinned  -Mcudalib=cufft -std=c++11" )
#   SET(CMAKE_EXEC_LINKER_FLAGS "${CMAKE_EXEC_LINKER_FLAGS} -O3  -Minfo=accel -ta=tesla,cuda9.0  -Mcudalib=cufft -std=c++11" )
  #  SET(CMAKE_EXEC_LINKER_FLAGS "${CMAKE_EXEC_LINKER_FLAGS} -g -Minfo=accel -ta=multicore  -Mcudalib=cufft -std=c++11" )
endif()


message(${RED} ${CMAKE_CXX_FLAGS} ${ColourReset} )
# Generate Object Files
add_executable(PittPack ${OBJ} ${HEADER_FILES}) 

#add_library()

# do not forget to link external libraries here 
#if(!GPU)
target_link_libraries(PittPack ${HDF5_LIBRARIES} ${FFTW_LIB} ${MPI_LIBRARIES} )
#endif
#target_link_libraries(amrGem )
#target_link_libraries(amrGem  )

#install(TARGETS amrGem DESTINATION ${dir/}bin)

############################################################################
#
#                    SECTION: ARCHIEVING
#
#	 1) it gets the date and appends that to the project name
#	 2) compresses the "/src" in a TGZ format
# 	 3) Puts the resulting file in the "/archieve" directory
#
############################################################################

MACRO (TODAY RESULT)
    IF (WIN32)
        EXECUTE_PROCESS(COMMAND "date" "/T" OUTPUT_VARIABLE ${RESULT})
         string(REGEX REPLACE "(..)/(..)/..(..).*" "\\2.\\1.\\3"
${RESULT} ${${RESULT}})
    ELSEIF(UNIX)
        EXECUTE_PROCESS(COMMAND "date" "+%d/%m/%Y" OUTPUT_VARIABLE ${RESULT})
        string(REGEX REPLACE "(..)/(..)/..(..).*" "\\2.\\1.\\3"
${RESULT} ${${RESULT}})
    ELSE (WIN32)
        MESSAGE(SEND_ERROR "date not implemented")
        SET(${RESULT} 000000)
    ENDIF (WIN32)
ENDMACRO (TODAY)

TODAY(DATE)
message(STATUS "Compilation date = ${DATE}")


# the following two lines specify the file type for binary and source
# generation separately

SET(CPACK_GENERATOR "TGZ")
SET(CPACK_SOURCE_GENERATOR "TGZ")

#INCLUDE(CPack)
# Do not rely on CPACK_xxx variables, they are not set till INCLUDE(CPack) is
# invoked, so there is no way to set them and as cpack to name the file as you
# like, I hard coded the project name and the date as a file_name to be used
# by cpack, this works fine 

#message("package name: " ${CPACK_PACKAGE_NAME})

set(CPACK_PACKAGE_FILE_NAME "PittPack-${DATE}")

set(CPACK_SOURCE_PACKAGE_FILE_NAME "PittPack-${DATE}")

# set the directory you would like to be compressed, For under development
# projects it is normally the src directory that includes *.cpp and *.h files

SET(CPACK_SOURCE_INSTALLED_DIRECTORIES "${dir}/src;/")

SET(CPACK_PACKAGE_VENDOR "PITT")

# If you like to add any other files to archieve, need to use install file
# command, now in addition to source I am compresing the CMakeList.txt along
# with itimply installing it to the source directory. since src is the primary
# directory to archieve, make install command will copy this to src directory
# though, use the  following line if you would like to include
# CMakeList.txt in your archieve otherwise comment it out, one way to 
# get rid of this is to use an alias to remove that file from src after
# packing

install(FILES ${dir}/CMakeLists.txt DESTINATION ${dir}/src)

message("file name:" ${CPACK_PACKAGE_FILE_NAME})

#install(TARGETS  amrGem ARCHIVE  DESTINATION archieve)

# This variable (CPACK_OUTPUT_FILE_PREFIX)  will put the result of amrGem.date.tar.gz into the folder
# named archieve, this is only available in  in CMake versions >= 2.8.3 

set(CPACK_OUTPUT_FILE_PREFIX "${dir}/archieve")

#message("PACKAGE Install Dir:" ${CPACK_PACKAGE_INSTALL_DIRECTORY})

INCLUDE(CPack)


#
#  Doxygen Integrated to Generate Documentation for the Code
#  Specify the INPUT parameter in Doxyfile as the src and include files
#  You need to tell doxygen where to look for files by giving the directory
#
# Do not forget to specify the input directory "../src/include" "../src" in
# Doxyfile

# st OFF to ON to enable documaneting
option(BUILD_DOCUMENTATION "Use Doxygen to create the HTML based API documentation" OFF)
if(BUILD_DOCUMENTATION)
  FIND_PACKAGE(Doxygen)

  if (NOT DOXYGEN_FOUND)
    message(FATAL_ERROR
      "Doxygen is needed to build the documentation. Please install it
correctly")
  endif()

  add_custom_target (doc ALL
                     ${DOXYGEN_EXECUTABLE}
${CMAKE_CURRENT_SOURCE_DIR}/doc/Doxyfile
                     WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/doc
                     )

endif()


