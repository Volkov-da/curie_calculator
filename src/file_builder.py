from pymatgen.core import Structure
from pymatgen.analysis.magnetism.analyzer import MagneticStructureEnumerator
from pymatgen.io.vasp.sets import MPRelaxSet, MPStaticSet
from shutil import copy
from tqdm import tqdm
import os

from pymatgen.io.vasp.outputs import Vasprun
import warnings
from time import sleep, gmtime, strftime
warnings.filterwarnings('ignore')


def create_job_script(out_path, ntasks=8, job_id=None):
    """
    Args:s
        out_path (str)   -   folder where job script will be created.
        job_id   (str)   -   preferable name of your job in squeue,
        and also the name of the folder for your vasp files.
        As a result in the folder with name 'job_id' will be created job script
    """
    if not job_id:
        job_id = os.path.basename(out_path)

    job_script_text = f"""#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks={ntasks}
#SBATCH --time=02:00:00
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
    out_data[6] = ' '.join(str(i)
                           for i in up_down_spin_counter(in_data)) + '\n'
    return out_data


def siman_POSCAR_writer(in_path: str, out_path: str) -> None:
    """
    Args:
        in_path  (str)  -   path to the POSCAR type file which needs to be made readable for siman
        out_path (str)  -   path where refactored version of this file will be written
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
    for folder in tqdm(afm_foldrs):
        tmp_out_path = os.path.join(
            out_path, 'POSCAR_' + folder.split("/")[-1])
        siman_POSCAR_writer(in_path=os.path.join(
            folder, 'POSCAR'), out_path=tmp_out_path)


def submit_all_jobs(input_folder: str) -> None:
    vasp_inputs_path = os.path.join(input_folder, 'vasp_inputs')
    initial_path = os.getcwd()
    for folder_name in os.listdir(vasp_inputs_path):
        os.chdir(initial_path)
        tmp_path = os.path.join(vasp_inputs_path, folder_name)
        os.chdir(tmp_path)
        os.system('sbatch jobscript.sh')
    os.chdir(initial_path)


LDAUJ_dict = {'Co': 0, 'Cr': 0, 'Fe': 0, 'Mn': 0, 'Mo': 0, 'Ni': 0, 'V': 0, 'W': 0,
              'Nb': 0, 'Sc': 0, 'Ru': 0, 'Rh': 0, 'Pd': 0, 'Cu': 0, 'Y': 0, 'Os': 0, 'Ti': 0, 'Zr': 0, 'Re': 0, 'Hf': 0, 'Pt': 0, 'La': 0}

LDAUU_dict = {'Co': 3.32, 'Cr': 3.7, 'Fe': 5.3, 'Mn': 3.9, 'Mo': 4.38, 'Ni': 6.2, 'V': 3.25, 'W': 6.2,
              'Nb': 1.45, 'Sc': 4.18, 'Ru': 4.29, 'Rh': 4.17, 'Pd': 2.96, 'Cu': 7.71, 'Y': 3.23, 'Os': 2.47, 'Ti': 5.89, 'Zr': 5.55,
              'Re': 1.28, 'Hf': 4.77, 'Pt': 2.95, 'La': 5.3}


LDAUL_dict = {'Co': 2, 'Cr': 2, 'Fe': 2, 'Mn': 2, 'Mo': 2, 'Ni': 2, 'V': 2, 'W': 2,
              'Nb': 2, 'Sc': 2, 'Ru': 2, 'Rh': 2, 'Pd': 2, 'Cu': 2, 'Y': 2, 'Os': 2,
              'Ti': 2, 'Zr': 2, 'Re': 2, 'Hf': 2, 'Pt': 2, 'La': 2}

relx_dict = {'ISMEAR': 0, 'SIGMA': 0.01, 'ISIF': 4, 'EDIFF': 1E-4, 'POTIM': 0.3,
             'EDIFFG': -0.01, 'SYMPREC': 1E-8, 'NCORE': 4, 'LCHARG': False, 'ICHARG': 2,
             'LDAU': True, 'LDAUJ': LDAUJ_dict, 'LDAUL': LDAUL_dict, 'LDAUU': LDAUU_dict, 'LWAVE': False,
             'LDAUPRINT': 1, 'LDAUTYPE': 2, 'LASPH': True, 'LMAXMIX': 4}

stat_dict = {'ISMEAR': -5, 'EDIFF': 1E-6, 'SYMPREC': 1E-8,  'NCORE': 4, 'ICHARG': 2,
             'LDAU': True, 'LDAUJ': LDAUJ_dict, 'LDAUL': LDAUL_dict, 'LDAUU': LDAUU_dict, 'NELM': 120, 'LVHAR': False,
             'LDAUPRINT': 1, 'LDAUTYPE': 2, 'LASPH': True, 'LMAXMIX': 4, 'LWAVE': False, 'LVTOT': False}


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
        vasp_relax_path (str)  -  path to the folder for static VASP run
        relax_dict      (dict) - dictionary with VASP INCAR keywords

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

    for i, magnetic_structure in tqdm(enumerate(enum_struct_list.ordered_structures)):

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


def static_changer(vasprun_path: str):
    """
    1. Replace INCAR relax for INCAR_stat
    2. Replace POSCAR with relaxed CONTCAR

    Namely prepare the folder for more accurate run of VASP
    for total energy estimation.
    """
    base_path = '/'.join(vasprun_path.split('/')[:-1])
    inc_path = os.path.join(base_path, 'INCAR')
    inc_stat_path = os.path.join(base_path, 'INCAR_stat')
    inc_relax_path = os.path.join(base_path, 'INCAR_relax')

    contcar_path = os.path.join(base_path, 'CONTCAR')
    poscar_path = os.path.join(base_path, 'POSCAR')

    log_relax = os.path.join(base_path, 'log_relax')
    log = os.path.join(base_path, 'log')

    out_relax = os.path.join(base_path, 'OUTCAR_relax')
    out = os.path.join(base_path, 'OUTCAR')

    copy(inc_path, inc_relax_path)  # INCAR -> INCAR_relax
    copy(inc_stat_path, inc_path)  # INCAR_stat -> INCAR
    copy(contcar_path, poscar_path)  # CONTCAR -> POSCAR
    copy(log, log_relax)  # log -> log_relax
    copy(out, out_relax)  # OUTCAR -> OUTCAR_relax


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
                    if vasprun.converged and vasprun.converged_ionic and vasprun.converged_electronic:
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
        sleep(20)
        if not vasprun_pathes:
            print(strftime("%H:%M:%S", gmtime()), 'All done!')
            break


def check_readiness(input_path: str, submit_path: str) -> None:
    submit_full_path = os.path.join(input_path, submit_path)
    pathes = sorted(os.listdir(submit_full_path))
    converged_list = []
    while len(converged_list) != len(pathes):
        tmp_time = strftime("%H:%M:%S", gmtime())
        print(f'{tmp_time} {submit_path.upper()} optimization')
        converged_list = []
        for folder_name in pathes:
            file_path = os.path.join(submit_full_path, folder_name, 'OSZICAR')
            if os.path.exists(file_path):
                with open(file_path) as f:
                    if 'E0' in f.readlines()[-1]:
                        converged_list += [folder_name]
                    else:
                        print(f'{folder_name} not yet converged')
            else:
                print(f'{folder_name} not yet written')
        print('\n')
        sleep(15)
    time = strftime("%H:%M:%S", gmtime())
    print(f'{time} {submit_path.upper()} optimization Finished')


def file_builder(input_path: str):
    assert os.path.exists(
        input_path), f'Input path: {input_path} does not exist!'
    assert os.path.exists(os.path.join(input_path, 'POSCAR')
                          ), f'Please specify POSCAR file in you input folder: {input_path}'

    get_VASP_inputs(input_path=input_path,
                    relx_dict=relx_dict,
                    static_dict=stat_dict)

    get_siman_inputs(input_path)
    copy(os.path.join(input_path, 'POSCAR'), os.path.join(
        input_path, 'siman_inputs', 'POSCAR_fm0'))
    print('All files written. Starting VASP calculations!')
    submit_all_jobs(input_path)
    sleep(7)
    check_readiness(input_path, submit_path='vasp_inputs')
    print('All AFM structures converged!. Starting J values estimation:')
    # vasprun_checker(input_path)


if __name__ == '__main__':
    input_path = os.getcwd()
    print(input_path)
    file_builder(input_path)
