#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=01:00:00
#SBATCH --job-name=tc_plot
#SBATCH --output=log
#SBATCH --error=err
#SBATCH -p lenovo

module load python/python36
python main_plotter.py