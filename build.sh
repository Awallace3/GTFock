#!/usr/bin/bash
rm -rf build
mkdir -p build
cd build
export CMAKE_PREFIX_PATH=$CONDA_PREFIX:$CMAKE_PREFIX_PATH
#export CC=icx
#export CXX=icpx
#export MPICC=mpiicx
#export MPICXX=mpiicpx
export CC=icc
export CXX=icpc
export MPICC=mpiicc
export MPICXX=mpiicpc

export AR=xiar rcs

cmake -DCMAKE_BUILD_TYPE=Debug -DCOMPILE_COMMANDS=ON ..
make
cd ..
