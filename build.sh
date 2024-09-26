#!/usr/bin/bash
# rm -rf build
# mkdir -p build
# export AR=xiar rcs

export CFLAGS=""

export SIMINT_LIBRARY_DIR=$PWD/../simint/build-avx512/install
export LIBCINT_LIBRARY_DIR=$PWD/../libcint/share/cmake/CInt
export GTMATRIX_LIBRARY_DIR=$PWD/../GTMatrix/share/cmake/GTMatrix
export ERD_LIBRARY_DIR=$PWD/../OptErd_Makefile/external/share/cmake/erd
export OED_LIBRARY_DIR=$PWD/../OptErd_Makefile/external/share/cmake/oed

export MPICC=$CONDA_PREFIX/bin/mpicc
export MPICXX=$CONDA_PREFIX/bin/mpicxx
export CC=$CONDA_PREFIX/bin/gcc
export CXX=$CONDA_PREFIX/bin/g++

export prefix_path="$SIMINT_LIBRARY_DIR"

# cmake -S. -Bbuild -DCMAKE_BUILD_TYPE=RelWithDebInfo -DCMAKE_EXPORT_COMPILE_COMMANDS=ON -DCInt_DIR=${LIBCINT_LIBRARY_DIR} -DGTMatrix_DIR=${GTMATRIX_LIBRARY_DIR} -Derd_DIR=$ERD_LIBRARY_DIR -Doed_DIR=$OED_LIBRARY_DIR -DCMAKE_C_COMPILER=$CC -DCMAKE_CXX_COMPILER=$CXX -G Ninja -DCMAKE_INSTALL_PREFIX=. -DCMAKE_PREFIX_PATH=${prefix_path}
ninja -C build 

# mpirun --mca btl ^openib -np 2 --oversubscribe build/pscf/pscf ./data/sto-3g.gbs ./data/water.xyz 2 1 1 4 10

# mpirun -np 2 ./build/pscf/pscf data/cc-pvdz/cc-pvdz.gbs data/alkane/alkane_62.xyz 2 1 1 2 10
mpirun -np 2 --oversubscribe build/pscf/pscf ./data/sto-3g.gbs ./data/water.xyz 2 1 1 4 10
