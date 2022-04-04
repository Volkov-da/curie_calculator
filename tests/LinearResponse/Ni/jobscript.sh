#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=04:00:00
#SBATCH --job-name=Lin_resp_Ni
#SBATCH --output=log
#SBATCH --error=err
#SBATCH -p lenovo

python ~/curie_calculator/src/linear_response.py 
