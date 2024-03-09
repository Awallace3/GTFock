#!/usr/bin/bash
rm -rf build
mkdir -p build
# export AR=xiar rcs

export SIMINT_LIBRARY_DIR=$PWD/../simint/build-avx512/install
export LIBCINT_LIBRARY_DIR=$PWD/../libcint/share/cmake/CInt
export GTMATRIX_LIBRARY_DIR=$PWD/../GTMatrix/share/cmake/GTMatrix

export MPICC=$CONDA_PREFIX/bin/mpicc
export MPICXX=$CONDA_PREFIX/bin/mpicxx
export CC=$CONDA_PREFIX/bin/gcc
export CXX=$CONDA_PREFIX/bin/g++

echo $MPICC
echo $MPICXX

# export prefix_path="$SIMINT_LIBRARY_DIR;${CONDA_PREFIX}/include"
export prefix_path="$SIMINT_LIBRARY_DIR"

cmake -S. -Bbuild -DCMAKE_EXPORT_COMPILE_COMMANDS=ON -DCInt_DIR=${LIBCINT_LIBRARY_DIR} -DGTMatrix_DIR=${GTMATRIX_LIBRARY_DIR} -DCMAKE_C_COMPILER=$CC -DCMAKE_CXX_COMPILER=$CXX -DCMAKE_MPICC_COMPILER=$MPICC -DCMAKE_MPICXX_COMPILER=$MPICXX -G Ninja -DCMAKE_INSTALL_PREFIX=. -DCMAKE_PREFIX_PATH=${prefix_path}
ninja -C build 
