#!/usr/bin/bash
# rm -rf build
mkdir -p build
# export AR=xiar rcs
export MPICC=$CONDA_PREFIX/bin/mpicc

export SIMINT_LIBRARY_DIR=$PWD/../simint/build-avx512/install
export LIBCINT_LIBRARY_DIR=$PWD/../libcint/share/cmake/CInt
export GTMATRIX_LIBRARY_DIR=$PWD/../GTMatrix/share/cmake/GTMatrix

export MPICXX=$CONDA_PREFIX/bin/mpicxx
cmake -S. -Bbuild -DCMAKE_EXPORT_COMPILE_COMMANDS=ON -DCInt_DIR=${LIBCINT_LIBRARY_DIR} -DCMAKE_PREFIX_PATH=${SIMINT_LIBRARY_DIR} -DGTMatrix_DIR=${GTMATRIX_LIBRARY_DIR} -DCMAKE_C_COMPILER=$CC -DCMAKE_CXX_COMPILER=$CXX -DCMAKE_MPICC_COMPILER=$MPICC -DCMAKE_MPICXX_COMPILER=$MPICXX -G Ninja -DCMAKE_INSTALL_PREFIX=.
ninja -C build 
