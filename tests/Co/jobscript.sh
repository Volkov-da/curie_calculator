#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=08:00:00
#SBATCH --job-name=Co
#SBATCH --output=log
#SBATCH --error=err
#SBATCH -p lenovo

python ../../src/stat_file_builder.py 