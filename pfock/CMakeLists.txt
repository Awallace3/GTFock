cmake_minimum_required(VERSION 3.10)
project(Pfock LANGUAGES C)

# Include the external make.in file here if needed, or ensure it's included from the parent CMakeLists.txt
# For example, you might have a CMakeLists.txt.in that you configure with configure_file()

# Set the C compiler to use MPI wrapper
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_C_COMPILER ${MPICC})

# Add compile options and definitions based on the configuration in make.in
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${CFLAGS}")

# Add preprocessor definitions
if(DYNAMIC EQUAL "1")
    add_definitions(-D__DYNAMIC__)
endif()

if(DO_SCF EQUAL "1")
    add_definitions(-D__SCF__)
endif()

# Include directories
# include_directories(
#     ${WORK_TOP}/libcint
#     ${WORK_TOP}/simint/build-avx512/install/include
#     ${WORK_TOP}/GTMatrix
#     ${WORK_TOP}/pscf
#     ${WORK_TOP}/pfock
# )

# Source files
file(GLOB SRC "*.c")

# Create a static library 'libpfock.a'
add_library(pfock STATIC ${SRC})

# Post-build command to copy the library and header file to a 'install' directory
add_custom_command(TARGET pfock POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E make_directory ../install
                   COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:pfock> ../install/
                   COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/pfock.h ../install/)


# Clean-up command (handled automatically by CMake)
set_property(DIRECTORY PROPERTY ADDITIONAL_MAKE_CLEAN_FILES "*.s;*.d;*~")


