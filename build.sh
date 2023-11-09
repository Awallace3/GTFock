#!/usr/bin/bash
rm -rf build
mkdir -p build
cd build
export CMAKE_PREFIX_PATH=$CONDA_PREFIX:$CMAKE_PREFIX_PATH
export CC=icx
export CXX=icpx
export AR=xiar rcs
export MPICC=mpiicx
export MPICXX=mpiicpx

cmake -DCMAKE_BUILD_TYPE=Debug -DCOMPILE_COMMANDS=ON ..
make
cd ..
