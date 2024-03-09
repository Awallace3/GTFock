cmake_minimum_required(VERSION 3.10)
project(Pscf LANGUAGES C)

# Assume the make.in file is included from the parent CMakeLists.txt

# Add compile options and definitions based on the configuration in make.in
# set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${CFLAGS}")
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -Wall -fopenmp -std=gnu99 -g -Wunknown-pragmas -Wunused-variable -I$ENV{CONDA_PREFIX}/include")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall -fopenmp -std=gnu99 -g -Wunknown-pragmas -Wunused-variable -I$ENV{CONDA_PREFIX}/include")

message(STATUS "PSCF CMAKE_C_FLAGS: ${CMAKE_C_FLAGS}")
message(STATUS "PSCF CMAKE_CXX_FLAGS: ${CMAKE_CXX_FLAGS}")

# Add preprocessor definitions
if(DYNAMIC EQUAL "1")
    add_definitions(-D__DYNAMIC__)
endif()

if(DO_SCF EQUAL "1")
    add_definitions(-D__SCF__)
endif()

# Source files for pscf
file(GLOB SRC "*.c")

# Libraries
set(BLAS_LIBS "-lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -lpthread -lm")
set(SCALAPACK_LIBS "-lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64")

add_executable(pscf ${SRC})

find_package(simint REQUIRED)
message(STATUS "SIMINT_FOUND: ${SIMINT_DIR}")
message(STATUS "SIMINT_LIBRARIES: ${simint_LIBRARIES}")
find_package(CInt REQUIRED)
find_package(GTMatrix REQUIRED)
find_package(MPI REQUIRED)

include_directories(
    ${MPI_C_INCLUDE_PATH} 
    ${simint_DIR}/include/simint
    ${CInt_INCLUDE_DIRS}
    ${GTMatrix_INCLUDE_DIRS}
    $ENV{CONDA_PREFIX}/include
    ${CMAKE_SOURCE_DIR}/aligned_malloc
)
# set(CMAKE_C_FLAGS " -O3 -Wall -qopenmp -std=gnu99 -fasm-blocks -g -xHost -Wunknown-pragmas -Wunused-variable  -diag-disable=10441")
message(STATUS "CMAKE_C_FLAGS: ${CMAKE_C_FLAGS}")
message(STATUS "simint_LIBS: ${simint_LIBRARIES}")
message(STATUS "CInt_LIBS: ${CInt_LIBRARIES}")
message(STATUS "MPI_C_LIBS: ${MPI_C_LIBRARIES}")
message(STATUS "CONDA_PREFIX: $ENV{CONDA_PREFIX}")

target_include_directories(pscf
    PUBLIC
    ${MPI_C_INCLUDE_PATH}
    ${simint_DIR}/include/simint
    ${CInt_INCLUDE_DIRS}
    ${GTMatrix_INCLUDE_DIRS}
    $ENV{CONDA_PREFIX}/include
    ${CMAKE_SOURCE_DIR}/aligned_malloc
)


target_link_libraries(pscf
    pfock
    aligned_malloc
    ${simint_LIBRARIES}
    ${CInt_LIBRARIES}
    ${GTMatrix_LIBRARIES}
    ${BLAS_LIBS}
    ${SCALAPACK_LIBS}
    ${MPI_C_LIBRARIES}
    $ENV{CONDA_PREFIX}/include
)



# Clean-up command (handled automatically by CMake)
set_property(DIRECTORY PROPERTY ADDITIONAL_MAKE_CLEAN_FILES "*.s;*.d;*~")
