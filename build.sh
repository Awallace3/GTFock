#!/usr/bin/bash
rm -rf build
mkdir -p build
cd build

export CMAKE_PREFIX_PATH=$CONDA_PREFIX:$CMAKE_PREFIX_PATH
export WORK_TOP=/theoryfs2/ds/amwalla3/gits

export CC=icx
export CXX=icpx
export MPICC=mpiicx
export MPICXX=mpiicpx
# export CC=icc
# export CXX=icpc
# export MPICC=mpiicc
# export MPICXX=mpiicpc

export AR=xiar rcs

cmake -DCMAKE_BUILD_TYPE=Debug -DCMAKE_EXPORT_COMPILE_COMMANDS=ON -DWORK_TOP=${WORK_TOP} ..
make
cd ..
