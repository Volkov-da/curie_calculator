import os
import numpy as np
from tqdm import tqdm
from shutil import copy, rmtree, move
from siman.calc_manage import smart_structure_read
from itertools import combinations
from scipy.constants import physical_constants
import matplotlib.pyplot as plt


PATH_TO_ENUMLIB = '../../enumlib'


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


def poscar_cleaner(in_data) -> list:
    """
    Args:
        in_data (list) - list of rows from any POSCAR type file

    Remove all unnecessary spaces from POSCAR files
    which makes them unreadable for other software (e.g siman)
    """
    out_data = []
    for st in in_data:
        if st.startswith('  '):
            st = st.replace('  ', '', 1)
            out_data.append(st)
        elif st.startswith(' '):
            st = st.replace(' ', '', 1)
            out_data.append(st)
        else:
            out_data.append(st)
    return out_data


def ratio_corrector(in_data: list) -> list:
    out_data = in_data.copy()
    stech_list = [int(i) for i in in_data[5].split()]
    corrected_stech = [sum(stech_list[:2]), *stech_list[2:]]
    ' '.join(str(i) for i in corrected_stech)
    out_data[5] = ' '.join(str(i) for i in corrected_stech) + '\n'
    out_data[6] = 'direct\n'
    return out_data


def atom_type_corrector(in_data: list, custom_atom_type=None) -> list:
    """
    Args:
        in_data           (list)  - list of lines from POSCAR type file.
        custom_atom_type  (str)   - string that you whant to write into the POSCAR
                                    as atom types e.g. "Fe O"
    """

    out_data = in_data.copy()
    if custom_atom_type:
        out_data[5] = custom_atom_type + '\n' + in_data[5]
    else:
        stech_str = in_data[0].split(' str #:')[0]
        right_atom_type = ''.join([i for i in stech_str if not i.isnumeric()])
        out_data[5] = right_atom_type + '\n' + in_data[5]
    return out_data


def poscar_pretiffier(in_path: str, out_path: str) -> None:
    """
    Args:
        in_path  (str)   - path to the POSCAR needed to be changed
        out_path (str)   - where to write changed POSCAR
    """
    with open(in_path) as in_f:
        in_data = in_f.readlines()

    out_data = atom_type_corrector(ratio_corrector(poscar_cleaner(in_data)))
    with open(out_path, 'w') as out_f:
        out_f.writelines(out_data)


def get_number_of_structures(enum_out='struct_enum.out') -> int:
    """
    Read file 'struct_enum.out' and
    returns number of generated supercells after
    the work of enum.x
    """
    if 'struct_enum.out' in os.listdir():
        with open(enum_out, "r") as file:
            lastline = (list(file)[-1])
            num_of_structures = int(lastline[:11])
        return num_of_structures
    else:
        print("\nERROR!\nWe need file 'struct_enum.out' to continue")


def direct_to_cart(in_path='CONTCAR', out_path='POSCAR_cart'):
    """

    This function transform CONTCAR file with direct
    coordinates into POSCAR with carthesian one.
    Take path to the CONCAR file as an input.
    Return file POSCAR_cart with carthesian coordiantes.

    """
    assert in_path in os.listdir(), f'{in_path} not here, nothing to transform!'
    with open(in_path, 'r') as contcar:
        contcar_text = contcar.readlines()

    head = contcar_text[: 7]
    lattice_param = np.loadtxt(contcar_text[2: 5])
    direct_coord = np.loadtxt(contcar_text[8:])
    cart_coord = direct_coord @ lattice_param

    np.savetxt('cart_coord.tmp', cart_coord)

    with open(out_path, 'w') as poscar_out:
        poscar_out.writelines(head)
        poscar_out.write('Carthesian\n')
        with open('cart_coord.tmp') as cart_coord:
            poscar_out.writelines(cart_coord.readlines())
        os.remove('cart_coord.tmp')


def run_enum(in_path: str) -> None:
    """
    Args:
        in_dir(str) - path to the folder where file "struct_enum.in" located
    Runs enumlib supercell generated based on prepare by user input file "struct_enum.in"
    """
    os.chdir(in_path)
    enum_exe_path = os.path.join(PATH_TO_ENUMLIB, 'src/enum.x')
    os.system(enum_exe_path)


def get_structures(path_to_enum=PATH_TO_ENUMLIB, num_of_structures=None) -> list:
    """
    This function read 'struct_enum.out'
    and generate POSCAR type files for all produced supercells.
    """
    if not num_of_structures:
        num_of_structures = get_number_of_structures()  # int number of generated supercells
    makeStrPath = os.path.join(PATH_TO_ENUMLIB, 'aux_src/makeStr.py')
    os.system(f'python {makeStrPath} 1 {num_of_structures}')
    print(f'Generated {num_of_structures} supercells')


def get_magmom_list(in_incar_data: list) -> list:
    magmom_line = [line for line in in_incar_data if 'MAGMOM' in line]
    magmom_list = [float(i) for i in magmom_line[0].split()[2:]]
    return magmom_list


def get_true_ratio(magmom_list: list, in_poscar_data: list) -> int:
    """
    Args:
        magmom_list    (list) - list of magnetic moments
        in_poscar_data (list) - list of string in correct POSCAR

    Return:
        true_ratio (int)
    """
    true_num_atoms = len(in_poscar_data[8:])
    current_num_atoms = len(magmom_list)
    true_ratio = int(true_num_atoms / current_num_atoms)
    return true_ratio


def magmom_lane_corrector(magmom_list: list, true_ratio: int) -> list:
    """
    Args:
        magmom_list (list)  - list of magnetic atoms from uncorrected INCAR file
        true_ratio  (int)   - ratio of uncorrected number atoms to the right number

    Returns:
        new_magmom_list (list) - list with corrected configuration for the afm cell

    Examples:
        magmom_lane_corrector([1, 1, 1, 0, 0], 2)
        >>> ([1, 1, 1, -1, -1, -1, 0, 0, 0, 0])

        magmom_lane_corrector([2, 0], 2)
        >>> ([2, -2, 0, 0])
    """

    magnetic_atoms_list = [i for i in magmom_list if i] * true_ratio
    noNmagnetic_atoms_list = [i for i in magmom_list if not i] * true_ratio
    middle_index = len(magnetic_atoms_list) // 2
    second_half = magnetic_atoms_list[middle_index:]
    magnetic_atoms_list[middle_index:] = [-i for i in second_half]
    new_magmom_list = magnetic_atoms_list + noNmagnetic_atoms_list
    return new_magmom_list


def incar_our_list_creator(in_incar_data: list, new_magmom_list: list) -> list:
    """
    Args:
        in_incar_data   (list)
        new_magmom_list (list)
    Returns:
        out_incar_data  (list) - list with lines for INCAR file, with correct MAGMOM line
    """
    for i, line in enumerate(in_incar_data):
        if 'MAGMOM' in line:
            magmom_line_index = i
    new_magmom_str = '    MAGMOM  =   ' + ' '.join([str(i) for i in new_magmom_list])
    out_incar_data = in_incar_data.copy()
    out_incar_data[magmom_line_index] = new_magmom_str
    return out_incar_data


def incar_pretiffier(in_path: str) -> None:
    """
    Args:
        in_path (str) -  path to the directory with INCAR file need to be corrected
        i.e. in_path = 'vasp_inputs/afm9/'
    """
    poscar_path = os.path.join(in_path, 'POSCAR')
    incar_path = os.path.join(in_path, 'INCAR')

    with open(poscar_path) as in_poscar:
        in_poscar_data = in_poscar.readlines()
    with open(incar_path) as in_incara:
        in_incar_data = in_incara.readlines()

    magmom_list = get_magmom_list(in_incar_data)
    true_ratio = get_true_ratio(magmom_list, in_poscar_data)
    new_magmom_list = magmom_lane_corrector(magmom_list, true_ratio)

    out_incar_data = incar_our_list_creator(in_incar_data, new_magmom_list)

    with open(incar_path, 'w') as out_f:
        out_f.writelines(out_incar_data)


def vasp_inputs_creator(num_of_structures: int, write_KPOINTS=False):
    """
    Args:
        num_of_structures  (int)  -   preferable number of structures you whant to study.
                                      20 is enough for majority of cases.
        write_KPOINTS      (bool) - write or not KPOINTS file from 'user_inputs' folders
                                    to all SP calculations for AFM structures.

    Creates folders ready for VASP calculations.
    Every folder contains one antiferromagnetic stucture that was generated by enumlib.
    """
    out_dir = 'vasp_inputs'
    input_directory = 'user_inputs'
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    current_path = os.getcwd()
    if not num_of_structures:
        num_of_structures = get_number_of_structures()
    for i in tqdm(range(num_of_structures)):
        tmp_path = os.path.join(current_path, out_dir, f'afm{i + 1}')
        os.makedirs(tmp_path, exist_ok=True)
        create_job_script(tmp_path, job_id=f'afm{i + 1}')
        poscar_pretiffier(in_path=f'vasp.{i + 1}', out_path=os.path.join(tmp_path, 'POSCAR'))
        copy(os.path.join(input_directory, 'INCAR'), os.path.join(tmp_path, 'INCAR'))
        incar_pretiffier(tmp_path)
        copy(os.path.join(input_directory, 'POTCAR'), os.path.join(tmp_path, 'POTCAR'))
        if write_KPOINTS:
            copy(os.path.join(input_directory, 'KPOINTS'), os.path.join(tmp_path, 'KPOINTS'))


def clean_all(input_folder: str) -> None:
    vasp_usless = [file for file in os.listdir(input_folder) if 'vasp.' in file]
    usless_files = ['debug_conc_check.out',
                    'debug_dvec_rots.out',
                    'debug_get_rotation_perms_lists.out',
                    'debug_site_restrictions.out',
                    'readcheck_enum.out',
                    'symops_enum_parent_lattice.out',
                    'VERSION.enum',
                    'struct_enum.out'
                    ] + vasp_usless
    for file in usless_files:
        try:
            os.remove(os.path.join(input_folder, file))
        except:
            continue

    rmtree(os.path.join(input_folder, 'vasp_inputs'), ignore_errors=True)
    rmtree(os.path.join(input_folder, 'siman_inputs'), ignore_errors=True)
    rmtree(os.path.join(input_folder, 'enum_out'), ignore_errors=True)


def afm_atoms_creator(in_data: list, custom_atom='Po') -> list:
    """
    Args:
        in_data (list) - list of rows from POSCAR type file.

    Add one type of "Fake" atom into the POSCAR structure.
    This allows it to be treated as an atoms with spin up and down respectively,
    thus estimate number of positive and negative contributions into the total energy.
    """
    out_data = in_data.copy()
    out_data[6] = 'direct\n'
    stech_str = in_data[0].split(' str #:')[0]
    right_atom_type = f'{custom_atom} ' + ''.join([i for i in stech_str if not i.isnumeric()])
    out_data[5] = right_atom_type + '\n' + in_data[5]
    return out_data


def siman_poscar_writer(in_path: str, out_path: str, custom_atom='Po') -> None:
    """
    Args:
        in_path  (str)  -   path to the POSCAR type file which needs to be made
                            readable for siman
        out_path (str)  -   path where refactored version of this file will be
                            written
    """

    with open(in_path) as in_f:
        in_data = in_f.readlines()

    out_data = afm_atoms_creator(poscar_cleaner(in_data), custom_atom=custom_atom)

    with open(out_path, 'w') as out_f:
        out_f.writelines(out_data)


def siman_inputs_creator(num_of_structures: int, out_dir='siman_inputs') -> None:
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    vasp_files = [file_name for file_name in os.listdir() if 'vasp.' in file_name]
    for i in range(1, num_of_structures + 1):
        siman_poscar_writer(f'vasp.{i}', os.path.join(out_dir, f'POSCAR_{i}'))


def check_magnetic_atoms(in_path: str) -> list:
    """
    Args:
        in_path (str) - path for the prepared to siman structure reader POSCAR type file.

    Returns:
        magnetic_atoms (list)  -    list of magnetic atoms in prepared
                                    structures with "fake" atoms
        e.g
        >>> magnetic_atoms
            ['Po', 'Eu']
    """
    with open(in_path, 'r') as in_f:
        in_data = in_f.readlines()
        magnetic_atoms = in_data[5].split()[: 2]
    return magnetic_atoms


def enum_out_collector(out_path='enum_out') -> None:
    if not os.path.exists(out_path):
        os.mkdir(out_path)
    files_list = [file for file in os.listdir() if os.path.isfile(file)]
    for file in files_list:
        if file == 'struct_enum.in':
            continue
        else:
            move(file, os.path.join(out_path, file))


def submit_all_jobs(input_folder: str) -> None:
    vasp_inputs_path = os.path.join(input_folder, 'vasp_inputs')
    initial_path = os.getcwd()

    for folder_name in os.listdir(vasp_inputs_path):
        if 'afm' in folder_name:
            os.chdir(initial_path)
            tmp_path = os.path.join(vasp_inputs_path, folder_name)
            os.chdir(tmp_path)
            os.system('sbatch jobscript.sh')
    os.chdir(initial_path)


def input_reader(path_to_input='./INPUT.txt') -> list:
    """
    Read INPUT.txt file and return input data
    Return:
        input_folder        (str)
        num_of_structures   (int)
        fake_magnetic_atoms (list)
        spin                (float)
    Return example:
        examples/Bi-Mn/
        14
        ['Po', 'Mn']
        2.5
    """
    with open(path_to_input) as in_f:
        in_data = in_f.readlines()
    inputs_dict = dict()
    for line in in_data:
        if 'fake_magnetic_atoms' in line:
            fake_magnetic_atoms = line.split('=')[1].split()
        elif 'input_folder' in line:
            input_folder = line.split('=')[1].split()[0]
        elif 'num_of_structures' in line:
            num_of_structures = int(line.split('=')[1])
        elif 'spin' in line:
            spin = float(line.split('=')[1])

    return input_folder, num_of_structures, fake_magnetic_atoms, spin
