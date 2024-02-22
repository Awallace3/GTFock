cmake_minimum_required(VERSION 3.10)
project(GTFock)

# Minimum version of CMake required
cmake_minimum_required(VERSION 3.10)

# Project name
project(MyProject)

# Set the compilers
# set(CMAKE_C_COMPILER "mpicc")
# set(CMAKE_CXX_COMPILER "mpicxx")
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Compiler options
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -Wall -qopenmp -std=gnu99 -fasm-blocks -g -xHost -Wunknown-pragmas -Wunused-variable  -diag-disable=10441")
# set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -Wall -fopenmp -std=gnu99 -g -march=native -Wunknown-pragmas -Wunused-variable")

# Libraries
set(BLAS_LIBS "-lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -lpthread -lm")
set(SCALAPACK_LIBS "-lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64")

find_package(simint REQUIRED)
find_package(CInt REQUIRED)
find_package(GTMatrix REQUIRED)
find_package(MPI REQUIRED)

# message(STATUS "CInt_INCLUDE_DIRS: ${CInt_INCLUDE_DIRS}")
# message(STATUS "CInt_LIBRARIES: ${CInt_LIBRARIES}")

# Include and link directories
include_directories(
    # ${WORK_TOP}/simint/build-avx512/install/includes
    # ${WORK_TOP}/libcint
    # ${WORK_TOP}/GTMatrix
    # ${WORK_TOP}/deps/include
    ${simint_DIR}/include/simint
    ${CInt_INCLUDE_DIRS}
    ${GTMatrix_INCLUDE_DIRS}
    ${CMAKE_CURRENT_SOURCE_DIR}/pscf
    ${CMAKE_CURRENT_SOURCE_DIR}/pfock
    ${CONDA_PREFIX}/include
)



#message(STATUS "MPI FOUND")
#message(STATUS "MPI_INCLUDE_PATH: ${MPI_INCLUDE_PATH}")
#include_directories(${MPI_C_INCLUDE_PATH})

# Add definitions for conditional compilation
if(DYNAMIC)
    add_definitions(-DDYNAMIC)
endif()
if(DO_SCF)
    add_definitions(-DDO_SCF)
endif()
if(SCF_TIMMING)
    add_definitions(-DSCF_TIMMING)
endif()
if(NOT SCF_OUT)
    add_definitions(-DSCF_OUT)
endif()

# Add subdirectories
add_subdirectory(pfock)
add_subdirectory(pscf)
