from pymatgen.io.vasp.sets import MPStaticSet
from stat_file_builder import create_job_script
from pymatgen.core import Structure
from shutil import copy as cp
from time import sleep
import os
from tabulate import tabulate
from pymatgen.io.vasp.outputs import Outcar

import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression
from time import sleep, gmtime, strftime
import matplotlib.pyplot as plt

import warnings
warnings.filterwarnings('ignore')


def get_incar(user_settings: dict) -> str:
    return '\n'.join([f'{k} = {v}' for k, v in user_settings.items()])


def write_incar(path: str, user_settings: dict):
    if not os.path.exists(path):
        os.mkdir(path)
    incar_path = os.path.join(path, 'INCAR')
    with open(incar_path, 'w+') as incar_file:
        incar_file.write(get_incar(user_settings))


def submit_job(path: str) -> None:
    initial_path = os.getcwd()
    os.chdir(path)
    os.system('sbatch jobscript.sh')
    os.chdir(initial_path)


def vasp_ready(oszicar_path: str) -> bool:
    if os.path.exists(oszicar_path):
        with open(oszicar_path) as f:
            text = f.readlines()
        if text and ('E0' in text[-1]):
            return True
        return False
    else:
        False


def write_calculation(path: str, structure, user_settings: dict, mode=None, mode_dict=None) -> None:
    if not os.path.exists(path):
        os.mkdir(path)
    settings = user_settings.copy()
    if mode == 'nscf':
        assert os.path.exists(os.path.join(path, 'CHGCAR')
                              ), f'ERROR: need CHGCAR in {path} for NSCF run'
        mode_dict.update({'ICHARG':  11})
        settings.update(mode_dict)
    elif mode == 'scf':
        settings.update(mode_dict)

    static_set = MPStaticSet(structure)
    static_set.get_vasp_input().write_input(path)
    write_incar(path, user_settings=settings)
    create_job_script(out_path=path, ntasks=8)
    submit_job(path)


def check_readiness(pathes: list) -> None:
    pathes = sorted(pathes)
    converged_list = []
    while len(converged_list) != len(pathes):
        print(f'\n{strftime("%H:%M:%S", gmtime())} Perturbations pipline')
        converged_list = []
        for folder_name in pathes:
            oszicar_path = os.path.join(folder_name, 'OSZICAR')
            if vasp_ready(oszicar_path):
                converged_list += [folder_name]
            else:
                print(f'\t{folder_name} not yet converged')
        print('\n')
        sleep(5)
    time = strftime("%H:%M:%S", gmtime())
    print(f'{time} Perturbations pipline finished!\n')


def get_d_el(outcar_path: str) -> float:
    try:
        outcar = Outcar(outcar_path)
        return outcar.charge[0]['d']
    except:
        return None


def get_results_df(scf_path_list: list, nscf_path_list: list) -> pd.DataFrame:
    scf_outcars = [os.path.join(path, 'OUTCAR') for path in scf_path_list]
    nscf_outcars = [os.path.join(path, 'OUTCAR') for path in nscf_path_list]

    # Number of electrons inf SCF and NSCF runs
    d_el_scf = [get_d_el(out) for out in scf_outcars]
    d_el_nscf = [get_d_el(out) for out in nscf_outcars]

    # distortions list gonna be the same in bouth cases SCF NSCF
    distortions = [float(scf_path.split('_')[-1])
                   for scf_path in scf_path_list]
    df_results = pd.DataFrame(
        {'SCF': d_el_scf, 'NSCF': d_el_nscf, 'v': distortions}).sort_values(by='v')
    return df_results


def plot_results(df_results: pd.DataFrame, input_path: str) -> float:
    df_results = df_results.dropna()
    # linear regression for slopes
    X_line = np.linspace(df_results['v'].min(),
                         df_results['v'].max()).reshape(-1, 1)
    X = df_results['v'].to_numpy().reshape(-1, 1)

    # estimate slope of scf run
    lr_scf = LinearRegression()
    y = df_results['SCF'].to_numpy().reshape(-1, 1)
    lr_scf.fit(X, y)
    y_scf = lr_scf.predict(X_line)
    slope_scf = lr_scf.coef_[0][0]

    # estimate slope of nscf run
    lr_nscf = LinearRegression()
    y = df_results['NSCF'].to_numpy().reshape(-1, 1)
    y_nscf = lr_nscf.fit(X, y).predict(X_line)
    slope_nscf = lr_nscf.coef_[0][0]

    # plot staff
    plt.figure(figsize=(8, 5), dpi=150)
    plt.scatter(df_results['v'], df_results['NSCF'],
                label=f'NSCF (slope {slope_nscf:.3f})', s=40)
    plt.plot(X_line, y_nscf, color='r')

    plt.scatter(df_results['v'], df_results['SCF'],
                label=f'SCF (slope {slope_scf:.3f})', marker='^', s=40)
    plt.plot(X_line, y_scf, color='b')
    U = 1 / slope_scf - 1 / slope_nscf
    plt.title(f'U={U:.3f}')

    plt.xticks(df_results['v'])

    plt.grid(alpha=.4)
    plt.tight_layout()
    plt.xlabel(r'$\alpha(eV)$')
    plt.ylabel('Number of d-electrons')
    plt.legend()

    output_path = os.path.join(input_path, 'output')
    if not os.path.exists(output_path):
        os.mkdir(output_path)

    plt.savefig(os.path.join(output_path, 'E_tot_plot.png'),
                bbox_inches='tight')
    return U


def linear_response(input_path: str):
    # WORK STARTS HERE
    structure = Structure.from_file(os.path.join(input_path, 'POSCAR'))
    # remove OSZICAR file is exists to track calculation readines
    calc_fold_path = os.path.join(input_path, 'CalcFold')
    if not os.path.exists(calc_fold_path):
        os.mkdir(calc_fold_path)
    groundstate_path = os.path.join(calc_fold_path, 'relax')
    gs_oszicar_path = os.path.join(groundstate_path, 'OSZICAR')
    gs_chgcar_path = os.path.join(groundstate_path, 'CHGCAR')

    if os.path.exists(gs_oszicar_path):
        os.remove(gs_oszicar_path)

    # groundstate calculation
    write_calculation(groundstate_path, structure, user_settings=STAT_DICT)
    while not vasp_ready(gs_oszicar_path):
        print('\tNot yet')
        sleep(5)
    print('Done')

    # necessary cause sometimes file CHCAR not rdy yet
    sleep(30)
    # perturbation stage
    scf_path_list, nscf_path_list = [], []

    distortions = np.arange(
        MIN_DISTORTION, MAX_DISTORTION+STEP_DISTORTION, STEP_DISTORTION).round(2)
    distortions = distortions[np.where(distortions != 0)]
    for v in distortions:
        MODE_DICT = {
            'LDAU':  True, 'LDAUTYPE':   3, 'LDAUL':   '2 -1 -1', 'LDAUPRINT':   2,
            'LDAUU':   f'{v} 0.00 0.00',
            'LDAUJ':   f'{v} 0.00 0.00'}

        # SCF calculations
        scf_path = os.path.join(calc_fold_path, f'scf_{v}')
        write_calculation(scf_path, structure,
                          user_settings=STAT_DICT, mode='scf', mode_dict=MODE_DICT)

        # NSCF calculations
        nscf_path = os.path.join(calc_fold_path, f'nscf_{v}')
        if not os.path.exists(nscf_path):
            os.mkdir(nscf_path)
        nscf_chgcar_path = os.path.join(nscf_path, 'CHGCAR')
        # copy charge density CHGCAR from the DFT groundstate calculation
        cp(gs_chgcar_path, nscf_chgcar_path)
        write_calculation(nscf_path, structure,
                          user_settings=STAT_DICT, mode='nscf', mode_dict=MODE_DICT)

        scf_path_list.append(scf_path)
        nscf_path_list.append(nscf_path)

    check_readiness(scf_path_list + nscf_path_list)

    print('Starting postprocessing.')
    # when calculations ready start postprocessing
    df_results = get_results_df(scf_path_list, nscf_path_list)
    print(tabulate(df_results, headers='keys', tablefmt='psql'))
    U = plot_results(df_results, input_path)
    print(f'\n\tEstimated U = {U:.3f}\n\nPostprocessing finished!')


if __name__ == '__main__':
    STAT_DICT = {
        'PREC':  'A', 'NCORE': 4,
        'EDIFF':  1E-6, 'ISMEAR':  0, 'SIGMA':  0.2,
        'ISPIN':  2, 'LORBIT':  11, 'LMAXMIX':  4, 'MAGMOM':  '5.0 5.0 5.0 5.0'
    }

    MIN_DISTORTION = -0.2
    MAX_DISTORTION = 0.2
    STEP_DISTORTION = 0.05

    linear_response(os.getcwd())
