from pymatgen.core import Structure
from pymatgen.analysis.magnetism.analyzer import MagneticStructureEnumerator, CollinearMagneticStructureAnalyzer
from pymatgen.io.vasp.inputs import Poscar, Incar, Kpoints
from pymatgen.io.vasp.sets import MPRelaxSet, MPStaticSet, MPSOCSet
import os

def create_job_script(out_path, job_id='JobName'):
    """
    Args:
        out_path (str)   -   folder where job script will be created.
        job_id   (str)   -   preferable name of your job in squeue,
        and also the name of the folder for your vasp files.
        As a result in the folder with name 'job_id' will be created job script
    """
    job_script_text = f"""#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --time=06:00:00
#SBATCH --job-name={job_id}
#SBATCH --output=log
#SBATCH --error=err
#SBATCH -p lenovo
module load mpi/impi-5.0.3 intel/mkl-11.2.3 vasp/vasp-5.4.4
mpirun vasp_std"""
    if not os.path.exists(out_path):
        os.mkdir(out_path)
    with open(f'{out_path}/jobscript.sh', 'w') as job:
        job.writelines(job_script_text)
        

def get_enum_structures(input_path : str, mag_prec: float=0.1, enum_prec : float=0.00001):
    init_structure = Structure.from_file(os.path.join(input_path, 'POSCAR'))
    enum_struct_list = MagneticStructureEnumerator(init_structure, transformation_kwargs={'symm_prec':mag_prec,'enum_precision_parameter':enum_prec})
    return enum_struct_list

def get_VASP_inputs(input_path : str, enum_struct_list : list, incar_dict : dict) -> None:
    for i, magnetic_structure in enumerate(enum_struct_list.ordered_structures):
        magnetic_type = enum_struct_list.ordered_structure_origins[i]
        str_id = magnetic_type + str(i)
        vasp_out_path = os.path.join(input_path, 'vasp_inputs', str_id)
        relax_set = MPStaticSet(structure = magnetic_structure, user_incar_settings = incar_dict ,reciprocal_density=1000,force_gamma=True)
        relax_set.get_vasp_input().write_input(vasp_out_path)
        create_job_script(vasp_out_path, job_id=str_id)

def afm_atom_creator(in_data: list, custom_atom='Po') -> list:
    """
    Args:
        in_data (list) - list of rows from POSCAR type file.

    Add one type of "Fake" atom into the POSCAR structure.
    This allows it to be treated as an atoms with spin up and down respectively,
    thus estimate number of positive and negative contributions into the total energy.
    """
    out_data = in_data.copy()
    out_data[5] = 'Po ' + out_data[5]
    return out_data

def up_down_spin_counter(in_data: list) -> list:
    spin_down = 0
    spin_up = 0
    no_spin = 0
    for row in in_data[8:]:
        if 'spin=-' in row:
            spin_down += 1
        elif 'spin=' in row:
            spin_up += 1
        else:
            no_spin += 1
    return [spin_up, spin_down, no_spin]

def spin_row_replacer(in_data: list) -> list:
    out_data = in_data.copy()
    out_data[6] = ' '.join(str(i) for i in up_down_spin_counter(in_data)) + '\n'
    return out_data


def siman_POSCAR_writer(in_path : str, out_path: str) -> None:
    """
    Args:
        in_path  (str)  -   path to the POSCAR type file which needs to be made
                            readable for siman
        out_path (str)  -   path where refactored version of this file will be
                            written
    """
    with open(in_path) as in_f:
        in_data = in_f.readlines()

    out_data = spin_row_replacer(afm_atom_creator(in_data))

    with open(out_path, 'w+') as out_f:
            out_f.writelines(out_data)

def get_siman_inputs(input_path: str):
    out_path = os.path.join(input_path, 'siman_inputs')
    if not os.path.exists(out_path):
        os.mkdir(out_path)
    vasp_inputs_path = os.path.join(input_path, 'vasp_inputs')
    afm_foldrs = [os.path.join(vasp_inputs_path, i)
                  for i in [i for i in os.listdir(vasp_inputs_path) if 'afm' in i]]
    for folder in afm_foldrs:
        tmp_out_path = os.path.join(out_path, 'POSCAR_' + folder.split("/")[-1])
        siman_POSCAR_writer(in_path=os.path.join(folder, 'POSCAR'), out_path=tmp_out_path)
        
def submit_all_jobs(input_folder: str) -> None:
    vasp_inputs_path = os.path.join(input_folder, 'vasp_inputs')
    initial_path = os.getcwd()
    for folder_name in os.listdir(vasp_inputs_path):
        os.chdir(initial_path)
        tmp_path = os.path.join(vasp_inputs_path, folder_name)
        os.chdir(tmp_path)
        os.system('sbatch jobscript.sh')
    os.chdir(initial_path)

def main_runner(input_path : str, incar_dict : dict):
    assert os.path.exists(input_path), f'Input path: {input_path} does not exist!'
    assert os.path.exists(os.path.join(input_path, 'POSCAR')), f'Please specify POSCAR file in you input folder: {input_path}'
    enum_struct_list = get_enum_structures(input_path)
    get_VASP_inputs(input_path, enum_struct_list, incar_dict)
    get_siman_inputs(input_path)
    submit_all_jobs(input_path)

input_path = '../examples/EuO_2/'


incar_dict = {'SYSTEM':"EuO",
'ISIF' : 2,
'NCORE':4,
'IBRION':2,
'EDIFF':1E-6,
'ENCUT':500,
'IBRION':2,
'LORBIT':11,
'LWAVE': False,
'LCHARG': False,
'LVHAR': False,
'LAECHG' : False, 
'ISMEAR':1,
'SIGMA':0.1}

main_runner(input_path, incar_dict)