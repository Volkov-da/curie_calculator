#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=12:00:00
#SBATCH --job-name=Fe2O3_monte_carlo
#SBATCH --output=log
#SBATCH --error=err
#SBATCH -p lenovo
../../../vampire/vampire-serial