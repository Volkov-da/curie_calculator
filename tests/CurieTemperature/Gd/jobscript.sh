#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=08:00:00
#SBATCH --job-name=Gd
#SBATCH --output=log
#SBATCH --error=err
#SBATCH -p lenovo

python ~/curie_calculator/src/stat_file_builder.py 