import os
import numpy as np
from tqdm import tqdm
from shutil import copy, rmtree, move
from siman.calc_manage import smart_structure_read
from itertools import combinations
from scipy.constants import physical_constants
import matplotlib.pyplot as plt


PATH_TO_ENUMLIB = '~/enumlib'


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


def poscar_cleaner(in_data):
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


def ratio_corrector(in_data):
    out_data = in_data.copy()
    stech_list = [int(i) for i in in_data[5].split()]
    corrected_stech = [sum(stech_list[:2]), *stech_list[2:]]
    ' '.join(str(i) for i in corrected_stech)
    out_data[5] = ' '.join(str(i) for i in corrected_stech) + '\n'
    out_data[6] = 'direct\n'
    return out_data


def atom_type_corrector(in_data, custom_atom_type=None):
    """
    Args:
        in_data(list) - list of lines from POSCAR type file.
        custom_atom_type(str) - string that you whant to write ito the POSCAR as atom types e.g. "Fe O"
    """

    out_data = in_data.copy()
    if custom_atom_type:
        out_data[5] = custom_atom_type + '\n' + in_data[5]
    else:
        stech_str = in_data[0].split(' str #:')[0]
        right_atom_type = ''.join([i for i in stech_str if not i.isnumeric()])
        out_data[5] = right_atom_type + '\n' + in_data[5]
    return out_data


def poscar_pretiffier(in_path, out_path):
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


def get_number_of_structures(enum_out='struct_enum.out'):
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


def run_enum(in_path):
    """
    Args:
        in_dir(str) - path to the folder where file "struct_enum.in" located
    Runs enumlib supercell generated based on prepare by user input file "struct_enum.in"
    """
    os.chdir(in_path)
    enum_exe_path = os.path.join(PATH_TO_ENUMLIB, 'src/enum.x')
    os.system(enum_exe_path)


def get_structures(path_to_enum=PATH_TO_ENUMLIB, num_of_structures=None):
    """
    This function read 'struct_enum.out'
    and generate POSCAR type files for all produced supercells.
    """
    if not num_of_structures:
        num_of_structures = get_number_of_structures()  # int number of generated supercells
    makeStrPath = os.path.join(PATH_TO_ENUMLIB, 'aux_src/makeStr.py')
    os.system(f'python {makeStrPath} 1 {num_of_structures}')
    print(f'Generated {num_of_structures} supercells')


def get_magmom_list(in_incar_data):
    magmom_line = [line for line in in_incar_data if 'MAGMOM' in line]
    magmom_list = [float(i) for i in magmom_line[0].split()[2:]]
    return magmom_list


def get_true_ratio(magmom_list, in_poscar_data):
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


def magmom_lane_corrector(magmom_list, true_ratio):
    """
    Args:
        magmom_list (list)  - list of magnetic atoms from uncorrected INCAR file
        true_ratio  (int)   - ratio of uncorrected number atoms to the right number

    Returns:
        new_magmom_list    (list) - list with corrected configuration for the afm cell

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


def incar_our_list_creator(in_incar_data, new_magmom_list):
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


def incar_pretiffier(in_path):
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


def vasp_inputs_creator(num_of_structures=None, out_dir='vasp_inputs', input_directory='user_inputs'):
    """
    Args:
        num_of_structures  (int)  -   preferable number of structures you whant to study.
                                      20 is enough for majority of cases.
        out_dir            (str)  -   name of folder you want to collec all the generated inputs
                                      for VASP.
        input_directory    (str)  -   folder with preferable inputs for VASP calculations.

    Creates folders ready for VASP calculations.
    Every folder contains one ferromagnetic stucture that was generated by enumlib.
    """
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    current_path = os.getcwd()
    copy(os.path.join(input_directory, 'all_jobs.sh'), out_dir)
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
        copy(os.path.join(input_directory, 'KPOINTS'), os.path.join(tmp_path, 'KPOINTS'))


def clean_all():
    vasp_usless = [file for file in os.listdir() if 'vasp.' in file]
    usless_files = ['debug_conc_check.out',
                    'debug_dvec_rots.out',
                    'debug_get_rotation_perms_lists.out',
                    'debug_site_restrictions.out',
                    'readcheck_enum.out',
                    'symops_enum_parent_lattice.out',
                    'VERSION.enum',
                    'struct_enum.out'
                    ] + vasp_usless
    for file in tqdm(usless_files):
        try:
            os.remove(file)
        except:
            continue
    rmtree('./vasp_inputs/', ignore_errors=True)
    rmtree('./siman_inputs/', ignore_errors=True)
    rmtree('./enum_out/', ignore_errors=True)


def afm_atoms_creator(in_data, custom_atom='Po'):
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


def siman_poscar_writer(in_path, out_path, custom_atom='Po'):
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


def siman_inputs_creator(num_of_structures, out_dir='siman_inputs'):
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    vasp_files = [file_name for file_name in os.listdir() if 'vasp.' in file_name]
    for i in range(1, num_of_structures):
        siman_poscar_writer(f'vasp.{i}', os.path.join(out_dir, f'POSCAR_{i}'))


def check_magnetic_atoms(in_path):
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


def enum_out_collector(out_path='enum_out'):
    if not os.path.exists(out_path):
        os.mkdir(out_path)
    files_list = [file for file in os.listdir() if os.path.isfile(file)]
    for file in files_list:
        if file == 'struct_enum.in':
            continue
        else:
            move(file, os.path.join(out_path, file))


def get_coefficients(poscars_to_check, up_to=2):
    """
    As a result this function will return coefficients,
    For E0, J1 and J2 for the Hisenber model in a form of matrix.
    Combaining obtained matrix with calculated energies (E_afm) of AFM configurations,
    One can calculate E0, J1, J2, solving a system of linear equations.

    Args:
        poscars_to_check - List of paths to POCARs files that should be checked.
        up_to - how many coordinations spheres to check, default check first two.

    Example:
        poscars_to_check = ['/EuO/POSCARS2siman/POSCAR_1',
                            '/EuO/POSCARS2siman/POSCAR_2',
                            '/EuO/POSCARS2siman/POSCAR_3']
        get_coefficients(poscars_to_check)
        >>> np.array([[0,  6],
                      [4, -6],
                      [-6,  0]]
    """
    nn_matrix = []
    for path_to_poscar in tqdm(poscars_to_check):
        # get a list with number of nearest
        nn_ls = list(count_nn(path_to_poscar, magnetic_atoms=FAKE_MAGNETIC_ATOMS).values())
        # neighbours for the structure
        nn_matrix.append(nn_ls)  # add list to he matrix
    nn_matrix = np.array(nn_matrix)
    return nn_matrix[..., :up_to]


def count_nn(path_to_poscar, magnetic_atoms):
    """
    calculated the number of nearest neighbors,
    to fit into the Heisenberg model.
    Get a path to POSCAR structure as an input,
    To avoid errors one should use prettified POSCAR,
    use poscar_prettifier() function first.

    Args:
        poscar_path(str) - path to the POSCAR file
        magnetic_atoms(list) - two type of atoms to be treated as a magnetic
                            with an opposite spins (up/down).
                            your POSCAR should contain this two types of atoms.
    Returns:
        dict(distance : number_of_neibours)
    """
    if not os.path.exists(path_to_poscar):
        print(f'File {path_to_poscar} does not exist!')
        return None
    st = smart_structure_read(path_to_poscar)
    st = st.replic([6, 6, 6])
    out = st.nn(i=1, n=50, silent=1)

    a = list(zip(out['el'][1:], out['dist'][1:]))
    # collect all the unique distances
    unique_dist = set(round(i[1], 3) for i in a if i[0] in magnetic_atoms)
    magnetic_dist_lst = [(el, round(dist, 3)) for el, dist in a if el in magnetic_atoms]

    dist_neighbNum = {_: 0 for _ in unique_dist}  # key- distane, value list of
    # neighbours at 1st 2nd 3d coordination spheres

    for dist in unique_dist:
        for el, distance in magnetic_dist_lst:
            if dist == distance:
                if el == magnetic_atoms[0]:
                    dist_neighbNum[dist] += 1
                elif el == magnetic_atoms[1]:
                    dist_neighbNum[dist] -= 1
    return dist_neighbNum


def get_coefficients(poscars_to_check, fake_magnetic_atoms, up_to=2):
    """
    As a result this function will return coefficients,
    For E0, J1 and J2 for the Hisenber model in a form of matrix.
    Combaining obtained matrix with calculated energies (E_afm) of AFM configurations,
    One can calculate E0, J1, J2, solving a system of linear equations.

    Args:
        poscars_to_check - List of paths to POCARs files that should be checked.
        up_to - how many coordinations spheres to check, default check first two.

    Example:
        poscars_to_check = ['/EuO/POSCARS2siman/POSCAR_1',
                            '/EuO/POSCARS2siman/POSCAR_2',
                            '/EuO/POSCARS2siman/POSCAR_3']
        get_coefficients(poscars_to_check)
        >>> np.array([[0,  6],
                      [4, -6],
                      [-6,  0]]
    """
    nn_matrix = []
    for path_to_poscar in tqdm(poscars_to_check):
        # get a list with number of nearest
        nn_ls = list(count_nn(path_to_poscar, magnetic_atoms=fake_magnetic_atoms).values())
        # neighbours for the structure
        nn_matrix.append(nn_ls)  # add list to he matrix
    nn_matrix = np.array(nn_matrix)
    return nn_matrix[..., :up_to]


def get_exchange_couplings(nn_matrix, energies_afm, spin):
    """
    Returns three float numbers: E0, J1, J2 respectively.
        E0 - is the part of total energy independent of the spin configuration
        J1, J2 - are the first and second nearest neighbor exchange constants

    Example:
        nn_matrix = np.array([[ 0  6],
                              [ 4 -6],
                              [-6  0]])

        energies_afm = np.array([-41.614048,
                                 -41.611114,
                                 -41.59929])

        get_exchange_coeff(nn_matrix, energies_afm)
        >>>-41.609258249999996 -0.00010548412698408951 -5.068518518515268e-05
    """
    nn_matrix = nn_matrix * spin * (spin + 1)
    nn_matrix = np.append(np.ones([len(nn_matrix), 1]), nn_matrix, 1)
    E0, J1, J2 = np.linalg.solve(nn_matrix, energies_afm)
    return E0, J1, J2


def calculate_Tc(J1, J2, z1, z2):
    k_B = physical_constants['Boltzmann constant in eV/K'][0]
    T_c = (J1 * z1 + J2 * z2) / (3 * k_B)
    return T_c


def total_num_neighbours(path_to_poscar, magnetic_atoms):
    """
    Return total number of magnetic_atoms in first,
    second and third coordination spheres z1, z2, z3 respectively.
    Since this number is constant for particular structure and will be the
    same for all generated supercells you can use whatever POSCAR file you want.

    Args:
        magnetic_atoms(list) - list of atoms which suppose to be magnetic in given structure.
        e.g [Fe, Co] or just [Fe] if there is only one type of magnetic atoms.

    Example:
        path_to_poscar = './vasp_inputs/afm3/POSCAR'
        total_num_neighbours(path_to_poscar, magnetic_atoms = ['Eu'])
        >>> (12, 6, 24)

    """
    if not os.path.exists(path_to_poscar):
        print(f'File {path_to_poscar} does not exist!')
    st = smart_structure_read(path_to_poscar)
    st = st.replic([6, 6, 6])
    out = st.nn(i=1, n=100, silent=1)

    a = list(zip(out['el'][1:], out['dist'][1:]))

    # collect all the unique distances
    unique_dist = set(round(i[1], 3) for i in a if i[0] in magnetic_atoms)
    magnetic_dist_lst = [round(dist, 3) for el, dist in a if el in magnetic_atoms]

    nn_num_lst = [magnetic_dist_lst.count(dist) for dist in unique_dist]

    z1 = nn_num_lst[0]
    z2 = nn_num_lst[1]
    z3 = nn_num_lst[2]
    return z1, z2, z3


def get_E_tot(in_path):
    """
    Args:
        in_path (str) - direct path to the log file from the VASP run
    Return:
        E_tot (float) - total calculated energies

    This function suppose to parse log type file from the VASP run
    and return the total energy of the structure (   1 F= -...)
    """
    with open(in_path) as in_f:
        text = in_f.readlines()
    lst_line = text[-1]
    E_tot = float(lst_line.split()[2])
    return E_tot


def get_all_energies(in_path, num_of_structures):
    """
    Args:
        in_path           (str)
        num_of_structures (int)
    Return:
        E_dict (dict) {id_of_structure : E_tot}

    Function parse all the folder with VASP calculations
    for AFM structures and returns the dictionary in format
    {id_of_structure : E_tot}
    """
    incar_path = os.path.join(in_path, 'vasp_inputs', f'afm1', 'INCAR')
    with open(incar_path) as in_data:
        incar_data = in_data.readlines()
    magmom_list = get_magmom_list(incar_data)
    E_dict = dict()
    for i in range(num_of_structures):
        tmp_path = os.path.join(in_path, 'vasp_inputs', f'afm{i + 1}')
        log_path = os.path.join(tmp_path, 'log')
        poscar_path = os.path.join(tmp_path, 'POSCAR')
        with open(poscar_path) as in_data:
            poscar_data = in_data.readlines()
        ratio = get_true_ratio(magmom_list, poscar_data)

        E = get_E_tot(log_path)
        E_dict[i + 1] = E / ratio
    return E_dict


def get_Tc_list(num_of_structures, nn_matrix, energies_afm, z1, z2, spin):
    all_combinations = list(combinations(range(0, num_of_structures), 3))
    Tc_list = []
    num_of_singular_matrix = 0

    for combination in all_combinations:
        matrix = nn_matrix[(combination), ...]
        energies = energies_afm[(combination), ...]
        try:
            E0, J1, J2 = get_exchange_couplings(matrix,
                                                energies,
                                                spin=spin)
            Tc = calculate_Tc(J1, J2, z1=z1, z2=z2)
            Tc_list.append(Tc)
        except:
            num_of_singular_matrix += 1
    return Tc_list, num_of_singular_matrix


def get_results(input_folder, num_of_structures, FAKE_MAGNETIC_ATOMS, SPIN):
    En_dict = get_all_energies(input_folder, num_of_structures)

    energies_afm = np.array(list(En_dict.values()), dtype='float')

    z1, z2, z3 = total_num_neighbours(path_to_poscar=os.path.join(input_folder, 'vasp_inputs/afm1/POSCAR'),
                                      magnetic_atoms=FAKE_MAGNETIC_ATOMS)

    siman_input_path = os.path.join(input_folder, 'siman_inputs')

    poscars_to_check = [os.path.join(
        siman_input_path, f'POSCAR_{i}') for i in range(1, num_of_structures + 1)]

    nn_matrix = get_coefficients(poscars_to_check, fake_magnetic_atoms=FAKE_MAGNETIC_ATOMS)

    Tc_list, num_of_singular_matrix = get_Tc_list(num_of_structures=num_of_structures,
                                                  nn_matrix=nn_matrix,
                                                  energies_afm=energies_afm,
                                                  z1=z1,
                                                  z2=z2,
                                                  spin=SPIN)
    out_dict = {'Tc_list': Tc_list,
                'num_of_singular_matrix': num_of_singular_matrix,
                'energies_afm': energies_afm,
                'nn_matrix': nn_matrix,
                'z1': z1,
                'z2': z2,
                'z3': z3}

    return out_dict


def plot_results(Tcs):
    min_Tc = min(Tcs)
    max_Tc = max(Tcs)
    points_num = len(Tcs)
    plt.figure(figsize=(8, 8), dpi=150)
    plt.scatter(range(points_num),
                Tcs,
                s=3,
                c='red')
    plt.grid(alpha=.4)
    plt.yticks(range(int(round(min_Tc - 100, -2)),
                     int(round(max_Tc + 100, -2)),
                     100), fontsize=7)
    plt.xticks([])
    plt.axhline(c='black')
    plt.xlim(0, points_num)
    plt.ylabel(r'$T_C, K$')
    plt.savefig('tc.pdf')
