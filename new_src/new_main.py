from pymatgen.core import Structure
from pymatgen.analysis.magnetism.analyzer import MagneticStructureEnumerator, CollinearMagneticStructureAnalyzer
from pymatgen.io.vasp.inputs import Poscar, Incar, Kpoints
from pymatgen.io.vasp.sets import MPRelaxSet, MPStaticSet, MPSOCSet
import os

from pymatgen.io.vasp.outputs import Vasprun
import numpy as np
import warnings
from time import sleep
warnings.filterwarnings('ignore')


def create_job_script(out_path, job_id=None):
    """
    Args:
        out_path (str)   -   folder where job script will be created.
        job_id   (str)   -   preferable name of your job in squeue,
        and also the name of the folder for your vasp files.
        As a result in the folder with name 'job_id' will be created job script
    """
    if not job_id:
        job_id = os.path.basename(out_path)

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


def siman_POSCAR_writer(in_path: str, out_path: str) -> None:
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


def write_static_set(structure, vasp_static_path: str, static_dict: dict) -> None:
    """
    Args:
        structure        (pymatgen.core.structure.Structure)
        vasp_static_path (str)  -  path to the folder for static VASP run
        static_dict      (dict) - dictionary with VASP INCAR keywords

    Write the following files into specified folder:
        INCAR_stat
        jobscript.sh
    """
    if not os.path.exists(vasp_static_path):
        os.mkdir(vasp_static_path)
    static_set = MPStaticSet(structure,
                             user_incar_settings=stat_dict,
                             reciprocal_density=300,
                             force_gamma=True)
    static_set.incar.write_file(os.path.join(vasp_static_path, 'INCAR_stat'))
    create_job_script(vasp_static_path)


def write_relax_set(structure, vasp_relax_path: str, relax_dict: dict) -> None:
    """
    Args:
        structure        (pymatgen.core.structure.Structure)
        vasp_static_path (str)  -  path to the folder for static VASP run
        static_dict      (dict) - dictionary with VASP INCAR keywords

    Write the following files into specified folder:
        INCAR
        POSCAR
        POTCAR
        KPOINTS
        jobscript.sh
    """
    if not os.path.exists(vasp_relax_path):
        os.mkdir(vasp_relax_path)
    relax_set = MPRelaxSet(structure=structure,
                           user_incar_settings=relx_dict,
                           user_kpoints_settings={'reciprocal_density': 300},
                           force_gamma=True)
    relax_set.get_vasp_input().write_input(vasp_relax_path)
    create_job_script(vasp_relax_path)


def get_VASP_inputs(input_path: str, relx_dict: dict, static_dict: dict) -> None:

    init_structure = Structure.from_file(os.path.join(input_path, 'POSCAR'))
    enum_struct_list = MagneticStructureEnumerator(init_structure,
                                                   transformation_kwargs={'symm_prec': 0.1,
                                                                          'enum_precision_parameter': 0.00001})

    if not os.path.exists(os.path.join(input_path, 'vasp_inputs')):
        os.mkdir(os.path.join(input_path, 'vasp_inputs'))

    for i, magnetic_structure in enumerate(enum_struct_list.ordered_structures):

        magnetic_type = enum_struct_list.ordered_structure_origins[i]

        str_id = magnetic_type + str(i)

        vasp_out_path = os.path.join(input_path, 'vasp_inputs', str_id)

        if not os.path.exists(vasp_out_path):
            os.mkdir(vasp_out_path)

        write_relax_set(structure=magnetic_structure,
                        vasp_relax_path=vasp_out_path,
                        relax_dict=relx_dict)

        write_static_set(structure=magnetic_structure,
                         vasp_static_path=vasp_out_path,
                         static_dict=stat_dict)


def vasprun_checker(input_path):
    vasp_inputs_path = os.path.join(input_path, 'vasp_inputs')
    vasprun_pathes = sorted([os.path.join(vasp_inputs_path, i, 'vasprun.xml')
                             for i in os.listdir(vasp_inputs_path)])
    tmp_vasprun = vasprun_pathes.copy()
    while 1:
        print(len(vasprun_pathes))
        for i, vasprun_path in enumerate(vasprun_pathes):
            print(i + 1, end=' ')
            if os.path.exists(vasprun_path):
                try:
                    vasprun = Vasprun(vasprun_path, parse_dos=False,
                                      parse_eigen=False, exception_on_bad_xml=False)
                    if vasprun.converged:
                        print(f'Converged! {vasprun_path}')
                        tmp_vasprun.remove(vasprun_path)
                    else:
                        print(f'Not converged! {vasprun_path}')
                        tmp_vasprun.remove(vasprun_path)
                except Exception:
                    print('Still running')
            else:
                print(f'{vasprun_path} not written yet!')

        vasprun_pathes = tmp_vasprun.copy()
        print('\n')
        sleep(5)
        if not vasprun_pathes:
            print('All done!')
            break


def main_runner(input_path: str):
    assert os.path.exists(input_path), f'Input path: {input_path} does not exist!'
    assert os.path.exists(os.path.join(input_path, 'POSCAR')
                          ), f'Please specify POSCAR file in you input folder: {input_path}'

    get_VASP_inputs(input_path=input_path,
                    relx_dict=relx_dict,
                    static_dict=stat_dict)

    get_siman_inputs(input_path)
    submit_all_jobs(input_path)
    vasprun_checker(input_path)


stat_dict = {'ISMEAR': -5,
             'EDIFF': 1E-6,
             'SYMPREC': 1E-8,
             'NCORE': 4,
             'NSIM': 4,
             'ICHARG': 2,
             'NELM': 120,
             'LVHAR': False,
             'LASPH': True,
             'LMAXMIX': 4,
             'LCHARG': False,
             'LWAVE': False,
             'LVTOT': False,
             'LAECHG': False}

relx_dict = {'ISIF': 4,
             'ISMEAR': 0,
             'SIGMA': 0.01,
             'EDIFF': 1E-4,
             'POTIM': 0.3,
             'EDIFFG': -0.01,
             'SYMPREC': 1E-8,
             'NCORE': 4,
             'NSIM': 4,
             'LCHARG': False,
             'ICHARG': 2,
             'LWAVE': False,
             'LASPH': True,
             'LMAXMIX': 4}

if __name__ == '__main__':
    input_path = '../examples/EuO_2/'
    main_runner(input_path)
