cmake_minimum_required(VERSION 3.10)
project(Pfock LANGUAGES C)

# Set the C compiler to use MPI wrapper
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Find MPI
find_package(MPI REQUIRED)

if(CMAKE_BUILD_TYPE STREQUAL "DEBUG" OR CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
  message(STATUS "DEBUG MODE")
  add_definitions(-DDEBUG)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -fopenmp -std=gnu99 -g -Wunknown-pragmas -Wunused-variable -I$ENV{CONDA_PREFIX}/include")
else()
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -Wall -fopenmp -std=gnu99 -g -Wunknown-pragmas -Wunused-variable -I$ENV{CONDA_PREFIX}/include -ftree-vectorize")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall -fopenmp -std=gnu99 -g -Wunknown-pragmas -Wunused-variable -I$ENV{CONDA_PREFIX}/include")
endif()

set(BLAS_LIBS "-lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -lpthread -lm")
set(SCALAPACK_LIBS "-lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64")

# Include directories
include_directories(
    ${MPI_C_INCLUDE_PATH}  # This line might now be redundant with target_include_directories
    ${simint_DIR}/include/simint
    ${CInt_INCLUDE_DIRS}
    ${GTMatrix_INCLUDE_DIRS}
    $ENV{CONDA_PREFIX}/include
    ${CMAKE_SOURCE_DIR}/aligned_malloc
)

# Add preprocessor definitions
if(DYNAMIC EQUAL "1")
    add_definitions(-D__DYNAMIC__)
endif()

if(DO_SCF EQUAL "1")
    add_definitions(-D__SCF__)
endif()

# Source files
file(GLOB SRC "*.c")

# Create a static library 'libpfock.a'
add_library(pfock STATIC ${SRC})

# Ensure MPI is properly linked
target_include_directories(pfock PUBLIC ${MPI_C_INCLUDE_PATH} {BLAS_LIBS} {SCALAPACK_LIBS})
target_link_libraries(pfock PUBLIC MPI::MPI_C ${BLAS_LIBS} ${SCALAPACK_LIBS})

# Post-build command to copy the library and header file to a 'install' directory
add_custom_command(TARGET pfock POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E make_directory ../install
                   COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:pfock> ../install/
                   COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/pfock.h ../install/
)

# Clean-up command (handled automatically by CMake)
set_property(DIRECTORY PROPERTY ADDITIONAL_MAKE_CLEAN_FILES "*.s;*.d;*~")
