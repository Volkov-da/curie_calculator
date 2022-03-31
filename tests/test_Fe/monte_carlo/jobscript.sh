#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=06:00:00
#SBATCH --job-name=test_Fe
#SBATCH --output=log
#SBATCH --error=err
#SBATCH -p lenovo
../vampire/vampire-serial