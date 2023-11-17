#!/usr/bin/bash
rm -rf build
mkdir -p build
cd build

export CMAKE_PREFIX_PATH=$CONDA_PREFIX:$CMAKE_PREFIX_PATH
export WORK_TOP=/theoryfs2/ds/amwalla3/projects/gtf_psi4

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
cmake -DCMAKE_BUILD_TYPE=Debug -DCMAKE_EXPORT_COMPILE_COMMANDS=ON -DWORK_TOP=${WORK_TOP} -DCInt_DIR=/theoryfs2/ds/amwalla3/projects/gtf_psi4/libcint -DGTMatrix_DIR=/theoryfs2/ds/amwalla3/projects/gtf_psi4/GTMatrix .. # -DCMAKE_C_COMPILER=icc -DCMAKE_CXX_COMPILER=icpc -DCMAKE_Fortran_COMPILER=ifort -DCMAKE_MPI_C_COMPILER=mpiicc -DCMAKE_MPI_CXX_COMPILER=mpiicpc
make
cd ..
