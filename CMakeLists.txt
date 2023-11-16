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
# For linking MKL and other libraries, you can use find_package if available, or set them manually
set(BLAS_LIBS "-lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -lpthread -lm")
set(SCALAPACK_LIBS "-lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64")
# set(SCALAPACK_LIBS "-L${MKLROOT}/lib/intel64 -lmkl_scalapack_ilp64 -Wl,--no-as-needed -lmkl_intel_ilp64 -lmkl_gnu_thread -lmkl_core -lmkl_blacs_openmpi_ilp64 -lgomp -lpthread -lm -ldl")



# Include and link directories
# include_directories(${MPI_INCDIR})
# link_directories(${MPI_LIBDIR})



# NOTE: UPDATE THIS TO YOUR OWN PATH
# set(WORK_TOP "/storage/home/hhive1/awallace43/gits")
message(STATUS "WORK_TOP: ${WORK_TOP}")

# Include and link directories
include_directories(
    ${WORK_TOP}/libcint
    ${WORK_TOP}/simint/build-avx512/install/include
    ${WORK_TOP}/GTMatrix
    ${CMAKE_CURRENT_SOURCE_DIR}/pscf
    ${CMAKE_CURRENT_SOURCE_DIR}/pfock
    ${CONDA_PREFIX}/include
)

find_package(MPI REQUIRED)

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
