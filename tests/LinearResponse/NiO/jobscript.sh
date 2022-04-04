#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=8:00:00
#SBATCH --job-name=LinResp_NiO
#SBATCH --output=log
#SBATCH --error=err
#SBATCH -p lenovo

. ~/curie_calculator/.venv/bin/activate
python ~/curie_calculator/src/linear_response.py 
