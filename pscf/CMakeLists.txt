cmake_minimum_required(VERSION 3.10)
project(Pscf LANGUAGES C)

# Assume the make.in file is included from the parent CMakeLists.txt

# Add compile options and definitions based on the configuration in make.in
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${CFLAGS}")
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Add preprocessor definitions
if(DYNAMIC EQUAL "1")
    add_definitions(-D__DYNAMIC__)
endif()

if(DO_SCF EQUAL "1")
    add_definitions(-D__SCF__)
endif()

# Source files for pscf
file(GLOB SRC "*.c")

# Create a static library 'libpscf.a'
# add_library(pscf STATIC ${SRC})
add_executable(pscf ${SRC})


target_link_libraries(pscf
    pfock
    # ${WORK_TOP}/libcint/libcint.a
    # ${WORK_TOP}/GTMatrix/build/libGTMatrix.a
    # ${WORK_TOP}/simint/build-avx512/install/lib/libsimint.a
    # ${WORK_TOP}/deps/lib/libcint.a
    # ${WORK_TOP}/deps/lib/libGTMatrix.a
    ${BLAS_LIBS}
    ${SCALAPACK_LIBS}
    ${MPI_C_LIBRARIES}
    )



# Clean-up command (handled automatically by CMake)
set_property(DIRECTORY PROPERTY ADDITIONAL_MAKE_CLEAN_FILES "*.s;*.d;*~")
