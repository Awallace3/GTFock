#!/usr/bin/bash
echo mpirun -np 4 build/pscf/pscf data/cc-pvdz/cc-pvdz.gbs data/alkane/alkane_62.xyz 2 2 2 4 10


mpirun -np 8 build/pscf/pscf data/cc-pvdz/cc-pvdz.gbs data/alkane/alkane_62.xyz 4 2 2 5 10
#mpirun -np 4 build/pscf/pscf data/cc-pvdz/cc-pvdz.gbs data/alkane/alkane_62.xyz 2 2 1 4 10
# mpirun -np 2 xterm -e 'gdb --args build/pscf/pscf data/cc-pvdz/cc-pvdz.gbs data/alkane/alkane_62.xyz 2 1 1 5 10'
