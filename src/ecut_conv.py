import os
import numpy as np
from time import sleep, gmtime, strftime
from pymatgen.core import Structure
from pymatgen.io.vasp.sets import MPStaticSet
from file_builder import create_job_script
from solver import energy_list_getter, find_good_structures
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
                         'LDAU': False,
                         'LVHAR': False,
                         'LCHARG': False,
                         'LAECHG': False,
                         'LASPH': False}
        ecut_path = os.path.join(input_path, 'encut', str(ecut))
        if not os.path.exists(ecut_path):
            os.mkdir(ecut_path)
        static_set = MPStaticSet(structure, user_incar_settings=user_settings)
        static_set.get_vasp_input().write_input(ecut_path)
        create_job_script(out_path=ecut_path, ntasks=24)


def en_per_atom_list(input_path: str) -> list:
    _, struct_list = find_good_structures(input_path, folder='encut')  # TODO: rewrite
    struct_list = sorted(struct_list)
    print(struct_list)
    ecut_opt_path = os.path.join(input_path, 'encut')
    folder_list = os.listdir(ecut_opt_path)
    initial_atoms_num = len(Structure.from_file(os.path.join(input_path, 'POSCAR')))
    en_tot_list = energy_list_getter(struct_list, initial_atoms_num)
    en_per_atom = en_tot_list / initial_atoms_num
    return en_per_atom


def get_ecut(en_per_atom: list, ecut_range) -> int:
    y = np.diff(en_per_atom * 1000)  # diff in meV
    x = ecut_range[1:]
    Ecut = x[np.argmin(abs(y))]
    return Ecut


def plot_encut(input_path: str, en_per_atom: list, ecut_range: list) -> None:
    x = ecut_range
    y = en_per_atom * 1000  # to meV
    Ecut = get_ecut(en_per_atom, ecut_range)
    i = list(ecut_range).index(Ecut)
    plt.figure(figsize=(12, 6), dpi=200)
    plt.scatter(x, y, c='b')
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
    for kpoints in kpoints_range:
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
        create_job_script(out_path=kpoints_path, ntasks=24)


def check_readiness(input_path: str, submit_path: str) -> None:
    submit_full_path = os.path.join(input_path, submit_path)
    initial_path = os.getcwd()
    pathes = sorted(os.listdir(submit_full_path))
    converged_list = []
    while len(converged_list) != len(pathes):
        print(strftime("%H:%M:%S", gmtime()))
        converged_list = []
        for folder_name in pathes:
            file_path = os.path.join(submit_full_path, folder_name, 'log')
            if os.path.exists(file_path):
                with open(file_path) as f:
                    if 'E0' in f.readlines()[-1]:
                        converged_list += [folder_name]
                    else:
                        print(f'{folder_name} not yet converged')
            else:
                print(f'{folder_name} not yet written')
        print()
        sleep(10)
    print(strftime("%H:%M:%S", gmtime()))
    print('All structures converged. Starting Encut estimation.')


if __name__ == '__main__':
    ecut_range = np.arange(400, 1000, 20)
    kpoints_range = np.arange(20, 150, 10)

    input_path = os.getcwd()
    print(input_path)
    get_ecut_files(input_path, ecut_range)
    submit_all_jobs(input_path=input_path, submit_path='encut')
    sleep(10)
    check_readiness(input_path, submit_path='encut')
    en_per_atom = en_per_atom_list(input_path)
    Ecut = get_ecut(en_per_atom, ecut_range)
    plot_encut(input_path, en_per_atom, ecut_range=ecut_range)
    print(f'Estimated value of Encut: {Ecut} eV\n')
    print('_' * 31)
    print(f'Starting KPOINTS optimization.')
    get_kpoints_files(input_path, Ecut, kpoints_range)
    submit_all_jobs(input_path=input_path, submit_path='kpoints')
    print()
    check_readiness(input_path, submit_path='kpoints')
    print(f'Kpoints optimization finished!')
