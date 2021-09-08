import os
import numpy as np
from tqdm import tqdm
from time import sleep, gmtime, strftime
from pymatgen.core import Structure
from pymatgen.io.vasp.sets import MPStaticSet
from pymatgen.io.vasp.outputs import Vasprun
from stat_file_builder import create_job_script
import matplotlib.pyplot as plt
plt.style.use('ggplot')


def submit_all_jobs(input_path: str, submit_path: str) -> None:
    submit_full_path = os.path.join(input_path, submit_path)
    initial_path = os.getcwd()
    for folder_name in os.listdir(submit_full_path):
        os.chdir(initial_path)
        tmp_path = os.path.join(submit_full_path, folder_name)
        os.chdir(tmp_path)
        os.system('sbatch jobscript.sh')
    os.chdir(initial_path)


def get_ecut_files(input_path: str, ecut_range: list) -> None:
    ecut_opt_path = os.path.join(input_path, 'encut')
    if not os.path.exists(ecut_opt_path):
        os.mkdir(ecut_opt_path)
    structure = Structure.from_file(os.path.join(input_path, 'POSCAR'))
    for ecut in ecut_range:
        user_settings = {'ENCUT': ecut, 'EDIFF': 1E-7, 'NCORE': 4,
                         'LDAU': False, 'LVHAR': False, 'LCHARG': False,
                         'LAECHG': False, 'LASPH': False}
        ecut_path = os.path.join(input_path, 'encut', str(ecut))
        if not os.path.exists(ecut_path):
            os.mkdir(ecut_path)
        static_set = MPStaticSet(structure, user_incar_settings=user_settings)
        static_set.get_vasp_input().write_input(ecut_path)
        create_job_script(out_path=ecut_path, ntasks=24)


def en_per_atom_list(input_path: str, mode: str) -> tuple:
    """

    Args:
        input_path (str): path to the folder with POSCAR file
        mode (str): type of calculation it might be 'encut' or 'kpoints'

    Returns:
        energy_arr: structures energies normalized by the number of atom
        params_range: range of optimized parameter (Encut, KPOINTS)
    """
    E_list = []
    calc_fold = os.path.join(input_path, mode)
    params_range = sorted([int(d) for d in os.listdir(calc_fold) if '.' not in d])
    for kpoint_density in tqdm(params_range):
        struct_folder = os.path.join(calc_fold, str(kpoint_density))
        vasprun_path = os.path.join(struct_folder, 'vasprun.xml')
        vasprun = Vasprun(vasprun_path, parse_dos=False, parse_eigen=False)
        E_tot = vasprun.final_energy
        E_list.append(E_tot)
    initial_atoms_num = len(Structure.from_file(os.path.join(input_path, 'POSCAR')))
    energy_arr = np.array(E_list) / initial_atoms_num
    return energy_arr, params_range


def get_ecut(en_per_atom: list, ecut_range) -> int:
    y = np.diff(en_per_atom * 1000)  # diff in meV
    x = ecut_range[1:]
    Ecut = x[np.argmin(abs(y))]
    return Ecut


def plot_encut(input_path: str, en_per_atom, ecut_range: list) -> None:
    x = ecut_range
    y = (en_per_atom - min(en_per_atom)) * 1000  # to meV
    Ecut = get_ecut(en_per_atom, ecut_range)
    i = list(ecut_range).index(Ecut)
    plt.figure(figsize=(12, 6), dpi=200)
    plt.plot(x, y, 'o-', c='r')
    plt.scatter(x[i], y[i], s=200, c='r')
    plt.ylabel('E/atom, meV')
    plt.xlabel('Encut, eV')
    plt.xticks(x, rotation=45, ha='right')
    plt.savefig(os.path.join(input_path, 'Encut.pdf'), bbox_inches='tight')


def write_kpoints(file_name, Rk):
    with open(file_name, 'w') as f:
        st = f"""Automatic mesh
0              ! number of k-points = 0 -> automatic generation scheme 
Auto           ! fully automatic
{Rk}           ! length (R_k)"""
        f.write(st)


def get_kpoints_files(input_path: str, Ecut: int, kpoints_range: list) -> None:
    kpoints_opt_path = os.path.join(input_path, 'kpoints')
    if not os.path.exists(kpoints_opt_path):
        os.mkdir(kpoints_opt_path)

    structure = Structure.from_file(os.path.join(input_path, 'POSCAR'))
    for kpoints in tqdm(kpoints_range):
        user_settings = {'ENCUT': Ecut,
                         'EDIFF': 1E-7, 'NCORE': 4,
                         'LDAU': False,
                         'LVHAR': False,
                         'LCHARG': False,
                         'LAECHG': False,
                         'LASPH': False}
        kpoints_path = os.path.join(input_path, 'kpoints', str(kpoints))
        file_name = os.path.join(kpoints_path, 'KPOINTS')
        if not os.path.exists(kpoints_path):
            os.mkdir(kpoints_path)
        static_set = MPStaticSet(structure, user_incar_settings=user_settings)
        static_set.get_vasp_input().write_input(kpoints_path)
        write_kpoints(file_name, Rk=kpoints)
        create_job_script(out_path=kpoints_path, ntasks=20)


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


def plot_kpoints(input_path: str, energy_arr, kpoints_range: list) -> None:
    y = (energy_arr - min(energy_arr)) * 1000
    x = kpoints_range
    plt.figure(figsize=(12, 6), dpi=200)
    plt.plot(x, y, 'o-', label='Sigma = 0.05', c='b')
    plt.ylabel('E/atom, meV')
    plt.xlabel(r'$R_K$')
    plt.xticks(x, rotation=45, ha='right')
    plt.legend()
    plt.savefig(os.path.join(input_path, 'Kgrid.pdf'), bbox_inches='tight')


def get_kpoint_density(energy_arr: list, kpoints_range) -> int:
    y = np.diff(energy_arr * 1000)  # diff in meV
    x = kpoints_range[1:]
    R_k = x[np.where(abs(y) < 0.5)[0][0]]
    return R_k


def encut_runner(input_path: str, ecut_range=np.arange(400, 900, 20)):
    print(f'Starting ENCUT optimization:')
    get_ecut_files(input_path, ecut_range)
    submit_all_jobs(input_path=input_path, submit_path='encut')
    check_readiness(input_path, submit_path='encut')
    encut_arr, encut_range = en_per_atom_list(input_path, mode='encut')
    Ecut = get_ecut(encut_arr, encut_range)
    plot_encut(input_path, encut_arr, encut_range)
    print(f'ENCUT optimization finished!')
    print('_' * 60 + '\n')
    return Ecut


def kpoints_runner(input_path: str, kpoints_range=np.arange(20, 100, 10)):
    print(f'Starting KPOINTS optimization:')
    get_kpoints_files(input_path, Ecut, kpoints_range)
    submit_all_jobs(input_path=input_path, submit_path='kpoints')
    print()
    check_readiness(input_path, submit_path='kpoints')
    print(f'KPOINTS optimization finished!')
    energy_arr, kpoints_range = en_per_atom_list(input_path, mode='kpoints')
    plot_kpoints(input_path, energy_arr, kpoints_range)
    R_k = get_kpoint_density(energy_arr, kpoints_range)
    print('_' * 60)
    return R_k


if __name__ == '__main__':
    input_path = os.getcwd()
    Ecut = encut_runner(input_path)
    print(f'Estimated value of Encut: {Ecut} eV\n')
    R_k = kpoints_runner(input_path)
    print(f'Estimated value of {R_k=}')
