from siman.calc_manage import smart_structure_read
import os
from tqdm import tqdm
from itertools import combinations
from scipy.constants import physical_constants
import numpy as np


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


def count_nn(path_to_poscar: str, magnetic_atoms: list) -> dict:
    """
    calculated the number of nearest neighbors,
    to fit into the Heisenberg model.
    Get a path to POSCAR structure as an input,
    To avoid errors one should use prettified POSCAR,
    use poscar_prettifier() function first.

    Args:
        poscar_path     (str)  - path to the POSCAR file
        magnetic_atoms  (list) - two type of atoms to be treated as a magnetic
                            with an opposite spins (up/down).
                            your POSCAR should contain this two types of atoms.
    Returns:
        dict{distance : number_of_neibours}
    """
    if not os.path.exists(path_to_poscar):
        print(f'File {path_to_poscar} does not exist!')
        return None
    st = smart_structure_read(path_to_poscar)
    st = st.replic([6, 6, 6])
    out = st.nn(i=1, n=500, silent=1)

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


def get_coefficients(poscars_to_check: list, magnetic_atoms: list, up_to=3) -> list:
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
        nn_ls = list(count_nn(path_to_poscar, magnetic_atoms=magnetic_atoms).values())
        # neighbours for the structure
        nn_matrix.append(nn_ls)  # add list to he matrix
    nn_matrix = np.array(nn_matrix)
    return nn_matrix[..., :up_to]


def get_exchange_couplings(nn_matrix: list, energies_afm: list, spin: float) -> list:
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
        >>>-41.609258249 0.000105484126 5.0685185185e-05
    """
    nn_matrix = nn_matrix * spin * (spin + 1)
    nn_matrix = np.append(np.ones([len(nn_matrix), 1]), nn_matrix, 1)
    E0, J1, J2, J3 = np.linalg.solve(nn_matrix, energies_afm)
    return E0, abs(J1), abs(J2), abs(J3)


def calculate_Tc(J1: float, J2: float, J3: float, z1: int, z2: int, z3: int) -> float:
    k_B = physical_constants['Boltzmann constant in eV/K'][0]
    T_c = 2 * (J1 * z1 + J2 * z2 + J3 * z3) / (3 * k_B)
    return T_c


def total_num_neighbours(path_to_poscar: str, magnetic_atoms: list) -> list:
    """
    Return total number of magnetic_atoms in first,
    second and third coordination spheres z1, z2, z3 respectively.
    Since this number is constant for particular structure and will be the
    same for all generated supercells you can use whatever POSCAR file you want.

    Args:
        path_to_poscar  (str)
        magnetic_atoms  (list)  - list of atoms which suppose to be magnetic in
        given structure.e.g [Fe, Co] or just [Fe] if there is only one type
        of magnetic atoms.

    Example:
        # BCC Iron case
        path_to_poscar = './vasp_inputs/afm1/POSCAR'
        total_num_neighbours(path_to_poscar, magnetic_atoms = ['Fe'])
        >>> [8, 6, 12, 24]

    """
    assert os.path.exists(path_to_poscar), f'File {path_to_poscar} does not exist!'
    num_neighb = [abs(i) for i in list(
        count_nn(path_to_poscar, magnetic_atoms=magnetic_atoms).values())]
    return num_neighb[:4]


def get_E_tot(in_path: str) -> float:
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


def get_all_energies(in_path: str, num_of_structures: int) -> dict:
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


def get_Tc_list(num_of_structures: int, nn_matrix, energies_afm: list, z1: int, z2: int, z3: int, spin: float):
    all_combinations = list(combinations(range(0, num_of_structures), 4))
    Tc_list = []
    Eg_list = []  # energy from the geometrly itself
    J1_list = []
    J2_list = []
    J3_list = []
    combination_list = []
    num_of_singular_matrix = 0

    for combination in all_combinations:
        matrix = nn_matrix[(combination), ...]
        energies = energies_afm[(combination), ...]
        try:
            Eg, J1, J2, J3 = get_exchange_couplings(matrix,
                                                    energies,
                                                    spin=spin)
            Tc = calculate_Tc(J1, J2, J3, z1=z1, z2=z2, z3=z3)
            Tc_list.append(Tc)
            Eg_list.append(Eg)
            J1_list.append(J1)
            J2_list.append(J2)
            J3_list.append(J3)
            combination_list.append(combination)
        except:
            num_of_singular_matrix += 1
    return Tc_list, Eg_list, J1_list, J2_list, J3_list, combination_list, num_of_singular_matrix


def get_results(input_folder: str, num_of_structures: int, magnetic_atoms: list, SPIN: float):
    En_dict = get_all_energies(input_folder, num_of_structures)

    energies_afm = np.array(list(En_dict.values()), dtype='float')

    z1, z2, z3, z4 = total_num_neighbours(path_to_poscar=os.path.join(input_folder, 'vasp_inputs/afm1/POSCAR'),
                                          magnetic_atoms=magnetic_atoms)

    siman_input_path = os.path.join(input_folder, 'siman_inputs')

    poscars_to_check = [os.path.join(
        siman_input_path, f'POSCAR_{i}') for i in range(1, num_of_structures + 1)]

    nn_matrix = get_coefficients(poscars_to_check, magnetic_atoms=magnetic_atoms)

    Tc_list, Eg_list, J1_list, J2_list, J3_list, combination_list, num_of_singular_matrix = get_Tc_list(num_of_structures=num_of_structures,
                                                                                                        nn_matrix=nn_matrix,
                                                                                                        energies_afm=energies_afm,
                                                                                                        z1=z1,
                                                                                                        z2=z2,
                                                                                                        z3=z3,
                                                                                                        spin=SPIN)
    out_dict = {'Tc_list': Tc_list,
                'num_of_singular_matrix': num_of_singular_matrix,
                'Eg_list': Eg_list,
                'J1_list': J1_list,
                'J2_list': J2_list,
                'J3_list': J3_list,
                'energies_afm': energies_afm,
                'nn_matrix': nn_matrix,
                'combination_list': combination_list,
                'z1': z1,
                'z2': z2,
                'z3': z3}

    return out_dict
