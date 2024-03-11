cmake_minimum_required(VERSION 3.10)
project(GTFock)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

find_package(MPI REQUIRED)
find_package(simint REQUIRED)
find_package(CInt REQUIRED)
find_package(GTMatrix REQUIRED)
find_package(erd REQUIRED)
find_package(oed REQUIRED)
find_package(BLAS REQUIRED)

# Compiler options
# set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -Wall -qopenmp -std=gnu99 -fasm-blocks -g -xHost -Wunknown-pragmas -Wunused-variable")
if(CMAKE_BUILD_TYPE STREQUAL "DEBUG" OR CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
  # Add debug-specific configurations
  message(STATUS "DEBUG MODE")
  add_definitions(-DDEBUG)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -fopenmp -std=gnu99 -g -Wunknown-pragmas -Wunused-variable -I$ENV{CONDA_PREFIX}/include")
else()
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -Wall -fopenmp -std=gnu99 -g -Wunknown-pragmas -Wunused-variable -I$ENV{CONDA_PREFIX}/include")
endif()
# set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D__SCF__")
message(STATUS "CMAKE_C_FLAGS: ${CMAKE_C_FLAGS}")

# Libraries
# set(BLAS_LIBS "-lmkl_gf_lp64 -lmkl_core -lmkl_gnu_thread -lpthread -lm")
set(BLAS_LIBS "-lblas -lmkl_gf_lp64 -lmkl_core -lmkl_gnu_thread -lm")
set(SCALAPACK_LIBS "-lmkl_scalapack_lp64 -lmkl_blacs_openmpi_lp64")

# Include and link directories
include_directories(
    ${MPI_C_INCLUDE_PATH}
    ${CMAKE_CURRENT_LIST_DIR}
    ${simint_DIR}/include/simint
    ${CInt_INCLUDE_DIRS}
    ${GTMatrix_INCLUDE_DIRS}
    $ENV{CONDA_PREFIX}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/aligned_malloc
    ${CMAKE_CURRENT_SOURCE_DIR}/pscf
    ${CMAKE_CURRENT_SOURCE_DIR}/pfock
    ${ERD_INCLUDE_DIRS}
    ${OED_INCLUDE_DIRS}
)



message(STATUS "MPI FOUND")
message(STATUS "MPI_INCLUDE_PATH: ${MPI_INCLUDE_PATH}")
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
add_subdirectory(aligned_malloc)
add_subdirectory(pfock)
add_subdirectory(pscf)
