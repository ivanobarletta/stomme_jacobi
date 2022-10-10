#!/bin/bash

#BSUB -q s_short
#BSUB -J stommel
#BSUB -n 1
#BSUB -o LOGS/stommel.%J.out
#BSUB -e LOGS/stommel.%J.err
#BSUB -P R000

module purge
module purge
module load intel19.5/19.5.281
module load impi19.5/19.5.281
module load impi19.5/petsc/3.7.5
module load intel19.5/netcdf/C_4.7.2-F_4.5.2_CXX_4.3.1

mpiexec.hydra -l ./stommel_jacobi

