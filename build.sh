#!/usr/bin/bash
# rm -rf build
mkdir -p build
cd build

export CPATH=/usr/include/x86_64-linux-gnu:$CPATH

export CMAKE_PREFIX_PATH=$CONDA_PREFIX:$CMAKE_PREFIX_PATH:/theoryfs2/ds/amwalla3/projects/gtf_psi4/simint/build-avx512/install
# export WORK_TOP=/theoryfs2/ds/amwalla3/projects/gtf_psi4

#export CC=icx
#export CXX=icpx
#export MPICC=mpiicx
#export MPICXX=mpiicpx
# export CC=icc
# export CXX=icpc
# export MPICC=mpiicc
# export MPICXX=mpiicpc

export AR=xiar rcs

# cmake -DCMAKE_BUILD_TYPE=Debug -DCMAKE_EXPORT_COMPILE_COMMANDS=ON -DWORK_TOP=${WORK_TOP} -DCInt_DIR=/theoryfs2/ds/amwalla3/projects/gtf_psi4/libcint -DGTMatrix_DIR=/theoryfs2/ds/amwalla3/projects/gtf_psi4/GTMatrix ..
# cmake -DCMAKE_BUILD_TYPE=Debug -DCMAKE_EXPORT_COMPILE_COMMANDS=ON -DCInt_DIR=/theoryfs2/ds/amwalla3/projects/gtf_psi4/libcint -DGTMatrix_DIR=/theoryfs2/ds/amwalla3/projects/gtf_psi4/GTMatrix .. # -DCMAKE_C_COMPILER=icc -DCMAKE_CXX_COMPILER=icpc -DCMAKE_Fortran_COMPILER=ifort -DCMAKE_MPI_C_COMPILER=mpiicc -DCMAKE_MPI_CXX_COMPILER=mpiicpc
# cmake -DCMAKE_BUILD_TYPE=Debug -DCMAKE_EXPORT_COMPILE_COMMANDS=ON -DCInt_DIR=/theoryfs2/ds/amwalla3/projects/gtf_psi4/libcint -DGTMatrix_DIR=/theoryfs2/ds/amwalla3/projects/gtf_psi4/GTMatrix .. # -DCMAKE_C_COMPILER=icc -DCMAKE_CXX_COMPILER=icpc -DCMAKE_Fortran_COMPILER=ifort -DCMAKE_MPI_C_COMPILER=mpiicc -DCMAKE_MPI_CXX_COMPILER=mpiicpc
# cmake -DCMAKE_BUILD_TYPE=Debug -DCMAKE_EXPORT_COMPILE_COMMANDS=ON -DCInt_DIR=/theoryfs2/ds/amwalla3/projects/gtf_psi4/libcint/share/cmake/CInt -DGTMatrix_DIR=/theoryfs2/ds/amwalla3/projects/gtf_psi4/GTMatrix/share/cmake/GTMatrix .. -DSIMINT_DIR=../simint/build-avx512/install
# cmake -DCMAKE_BUILD_TYPE=Debug -DCMAKE_EXPORT_COMPILE_COMMANDS=ON -DCInt_DIR=/theoryfs2/ds/amwalla3/projects/gtf_psi4/libcint/share/cmake/CInt -DGTMatrix_DIR=/theoryfs2/ds/amwalla3/projects/gtf_psi4/GTMatrix/share/cmake/GTMatrix .. -DSIMINT_DIR=../simint/build-avx512/install
# cmake -DCMAKE_BUILD_TYPE=Debug -DCMAKE_EXPORT_COMPILE_COMMANDS=ON -DCInt_DIR=/theoryfs2/ds/amwalla3/projects/gtf_psi4/libcint/share/cmake/CInt -DGTMatrix_DIR=/theoryfs2/ds/amwalla3/projects/gtf_psi4/GTMatrix/share/cmake/GTMatrix .. -Dsimint_DIR=../simint/build-avx512/install -DCMAKE_C_COMPILER=icc -DCMAKE_CXX_COMPILER=icpc -DCMAKE_Fortran_COMPILER=ifort
#cmake -DCMAKE_BUILD_TYPE=Debug -DCMAKE_EXPORT_COMPILE_COMMANDS=ON -DCInt_DIR=/theoryfs2/ds/amwalla3/projects/gtf_psi4/libcint/share/cmake/CInt -DGTMatrix_DIR=/theoryfs2/ds/amwalla3/projects/gtf_psi4/GTMatrix/share/cmake/GTMatrix .. -DCMAKE_C_COMPILER=icx -DCMAKE_CXX_COMPILER=icpx -DCMAKE_Fortran_COMPILER=ifort -DCMAKE_PREFIX_PATH=${CMAKE_PREFIX_PATH}
#cmake -DCMAKE_BUILD_TYPE=Debug -DCMAKE_EXPORT_COMPILE_COMMANDS=ON -DCInt_DIR=/theoryfs2/ds/amwalla3/projects/gtf_psi4/libcint/share/cmake/CInt -DGTMatrix_DIR=/theoryfs2/ds/amwalla3/projects/gtf_psi4/GTMatrix/share/cmake/GTMatrix ..  -DCMAKE_Fortran_COMPILER=ifort -DCMAKE_PREFIX_PATH=${CMAKE_PREFIX_PATH}
cmake -DCMAKE_BUILD_TYPE=Debug -DCMAKE_EXPORT_COMPILE_COMMANDS=ON -DCInt_DIR=/theoryfs2/ds/amwalla3/projects/gtf_psi4/libcint/share/cmake/CInt -DGTMatrix_DIR=/theoryfs2/ds/amwalla3/projects/gtf_psi4/GTMatrix/share/cmake/GTMatrix .. -DCMAKE_C_COMPILER=icc -DCMAKE_CXX_COMPILER=icpc -DCMAKE_Fortran_COMPILER=ifort -DCMAKE_PREFIX_PATH=${CMAKE_PREFIX_PATH}
make
cd ..
