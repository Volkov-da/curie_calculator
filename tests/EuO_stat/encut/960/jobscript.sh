#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --time=06:00:00
#SBATCH --job-name=960
#SBATCH --output=log
#SBATCH --error=err
#SBATCH -p lenovo
module load mpi/impi-5.0.3 intel/mkl-11.2.3 vasp/vasp-5.4.4
mpirun vasp_std