#!/usr/bin/bash
# mpirun -np 8 build/pscf/pscf data/cc-pvdz/cc-pvdz.gbs data/alkane/alkane_62.xyz 4 2 2 5 10
# srun -n 8 build/pscf/pscf data/cc-pvdz/cc-pvdz.gbs data/alkane/alkane_62.xyz 4 2 2 5 10
# export ORTE_BASE_HELP_AGGREGATE=0
# srun -N 1 -n 2  \
# mpirun --mca btl ^openib -np 2 build/pscf/pscf data/cc-pvdz/cc-pvdz.gbs data/alkane/alkane_62.xyz 2 1 1 5 10
# mpirun --mca btl ^openib -np 8 build/pscf/pscf data/cc-pvdz/cc-pvdz.gbs data/alkane/alkane_62.xyz 4 2 2 5 10

# ALKANE example
# mpirun --mca btl ^openib -np 8 --oversubscribe build/pscf/pscf data/cc-pvdz/cc-pvdz.gbs data/alkane/alkane_62.xyz 4 2 2 5 10

# WATER EXAMPLE
# mpirun --mca btl ^openib -np 2 --oversubscribe build/pscf/pscf ./data/aug-cc-pvtz.gbs ./data/water.xyz 2 1 1 4 10
# mpirun --mca btl ^openib -np 2 --oversubscribe build/pscf/pscf ./data/sto-3g.gbs ./data/water.xyz 2 1 1 4 10
mpirun  -np 2 build/pscf/pscf ./data/sto-3g.gbs ./data/water.xyz 2 1 1 4 10
