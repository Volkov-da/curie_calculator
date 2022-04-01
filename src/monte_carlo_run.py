import os

VAMPIRE_PATH = '../../../vampire/vampire-serial'

def input_monte_carlo(input_path: str, MAX_T: int, UNIT_CELL: float) -> None:
    out_path = os.path.join(input_path, 'monte_carlo')
    if not os.path.exists(out_path):
        os.mkdir(out_path)
    input_text = f"""#------------------------------------------
# Creation attributes:
#------------------------------------------
create:crystal-structure=sc
#------------------------------------------
# System Dimensions:
#------------------------------------------
dimensions:unit-cell-size = {UNIT_CELL} !A
dimensions:system-size-x = 7.7 !nm
dimensions:system-size-y = 7.7 !nm
dimensions:system-size-z = 7.7 !nm
#------------------------------------------
# Material Files:
#------------------------------------------
material:file=structure.mat
#------------------------------------------
# Simulation attributes:
#------------------------------------------
sim:minimum-temperature=0.0
sim:maximum-temperature={MAX_T}
sim:temperature-increment=20
sim:time-steps-increment=100
sim:equilibration-time-steps=10000
sim:loop-time-steps=30000
sim:time-step=1.0E-16
sim:applied-field-strength=0.0 !T
sim:applied-field-unit-vector=1,0,0
#------------------------------------------
# Program and integrator details
#------------------------------------------
sim:program=curie-temperature
sim:integrator=llg-heun
#------------------------------------------
# data output
#------------------------------------------
output:real-time
output:temperature
output:magnetisation
output:output-rate = 50
screen:time-steps
screen:temperature
screen:magnetisation-length
screen:mean-magnetisation-length"""
    file_out_path = os.path.join(out_path, 'input')
    with open(file_out_path, 'w+') as out_f:
        out_f.writelines(input_text)


def structure_monte_carlo(input_path: str, J_val: float, magmom: float):
    out_path = os.path.join(input_path, 'monte_carlo')
    if not os.path.exists(out_path):
        os.mkdir(out_path)

    input_text = f"""#---------------------------------------------------
# Number of Materials
#---------------------------------------------------
material:num-materials=1
# material[1]:material-name=Co
material[1]:damping-constant=1
material[1]:exchange-matrix[1]={J_val}
material[1]:atomic-spin-moment={magmom} !muB
material[1]:minimum-height=0.0
material[1]:maximum-height=1
material[1]:initial-spin-direction=0,0,1"""
    file_out_path = os.path.join(out_path, 'structure.mat')
    with open(file_out_path, 'w+') as out_f:
        out_f.writelines(input_text)


def job_monte_carlo(input_path: str, vampire_path: str, job_id=None) -> None:
    if not job_id:
        job_id = os.path.basename(input_path)

    out_path = os.path.join(input_path, 'monte_carlo')
    job_script_text = f"""#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=06:00:00
#SBATCH --job-name={job_id}
#SBATCH --output=log
#SBATCH --error=err
#SBATCH -p lenovo
{vampire_path}"""
    if not os.path.exists(input_path):
        os.mkdir(input_path)
    file_out_path = os.path.join(out_path, 'jobscript.sh')
    with open(file_out_path, 'w') as job:
        job.writelines(job_script_text)


def submit_monte_carlo(input_path: str) -> None:
    initial_path = os.getcwd()
    tmp_path = os.path.join(input_path, 'monte_carlo')
    os.chdir(tmp_path)
    os.system('sbatch jobscript.sh')
    os.chdir(initial_path)


def run_monte_carlo(input_path: str, MAX_T) -> None:
    input_monte_carlo(input_path, MAX_T, UNIT_CELL)
    structure_monte_carlo(input_path, J_val, magmom)
    job_monte_carlo(input_path, vampire_path=VAMPIRE_PATH)
    submit_monte_carlo(input_path)


UNIT_CELL = 3.54
MAX_T = 1400  # maximum temperature in MC simulation
J_val = 11.2e-21
magmom = 1.72

if __name__ == '__main__':
    run_monte_carlo(os.getcwd(), MAX_T)
