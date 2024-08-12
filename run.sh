#/usr/bin/bash
# mpirun -np 8 build/pscf/pscf data/cc-pvdz/cc-pvdz.gbs data/alkane/alkane_62.xyz 4 2 2 5 10
# srun -n 8 build/pscf/pscf data/cc-pvdz/cc-pvdz.gbs data/alkane/alkane_62.xyz 4 2 2 5 10
# export ORTE_BASE_HELP_AGGREGATE=0
# srun -N 1 -n 2  \
mpirun -np 2  \
    build/pscf/pscf data/cc-pvdz/cc-pvdz.gbs data/alkane/alkane_62.xyz 2 1 1 5 10
    # --mca orte_base_help_aggregate 0 \
# --mca btl_base_verbose 100 --mca btl_openib_verbose 100 \
