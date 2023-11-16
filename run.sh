#!/usr/bin/bash
conda deactivate
. /theoryfs2/ds/amwalla3/intel/oneapi/setvars.sh
export UCX_TLS=ud,sm,self
mpirun -np 8 pscf/scf data/cc-pvdz/cc-pvdz.gbs data/alkane/alkane_62.xyz 4 2 2 5 10
