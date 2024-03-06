cmake_minimum_required(VERSION 3.10)
project(aligned_malloc LANGUAGES C)

# Set the C compiler to use MPI wrapper
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

# Source files
file(GLOB SRC "*.c")
message(STATUS "Creating aligned_malloc")

# Create a static library 'libaligned_malloc.a'
add_library(aligned_malloc STATIC ${SRC})

# Post-build command to copy the library and header file to a 'install' directory
# Specify the library's include directories
target_include_directories(aligned_malloc PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
add_custom_command(TARGET aligned_malloc POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E make_directory ../install
                   COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:aligned_malloc> ../install/
                   COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/aligned_malloc.h ../install/)

# Clean-up command (handled automatically by CMake)
set_property(DIRECTORY PROPERTY ADDITIONAL_MAKE_CLEAN_FILES "*.s;*.d;*~")
# print out where install directory is
message(STATUS "Install directory: ${CMAKE_CURRENT_SOURCE_DIR}/../install")

