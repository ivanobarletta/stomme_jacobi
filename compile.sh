#!/bin/bash

#BSUB -q s_short
#BSUB -J mk
#BSUB -n 1
#BSUB -o logcompile.out
#BSUB -e logcompile.err
#BSUB -P R000

module purge
module purge
module load intel19.5/19.5.281
module load impi19.5/19.5.281
module load impi19.5/petsc/3.7.5

rm -f logcompile.*

mpiifort -I/zeus/opt/intel19.5/netcdf/C_4.7.2-F_4.5.2_CXX_4.3.1/include  -I/zeus/opt/impi19.5/petsc/3.7.5/include -I/zeus/opt/impi19.5/petsc/3.7.5/include/petsc/finclude -L/zeus/opt/impi19.5/petsc/3.7.5/lib -L/zeus/opt/intel19.5/netcdf/C_4.7.2-F_4.5.2_CXX_4.3.1/lib -lnetcdff -lpetsc main.f90 -o stommel_jacobi
