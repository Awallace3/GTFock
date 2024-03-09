#/usr/bin/bash
mpirun -np 8 build/pscf/pscf data/cc-pvdz/cc-pvdz.gbs data/alkane/alkane_62.xyz 4 2 2 5 10
