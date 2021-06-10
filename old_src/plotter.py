import os
import numpy as np
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")

"""
List of functions:
    1. plot_TCs(TCs:list, input_folder: str)

    2. plot_E_tot(energies_afm : list, nn_matrix : list, input_folder: str)

    3. plot_E_geom(Eg_list: list, combination_list : list, input_folder : str)
"""
DPI = 150
SIZE = (10, 8)


def plot_TCs(TCs: list, input_folder: str):
    TCs = np.array(TCs)
    text = f'mean: {TCs.mean():.2f} K\nmax  : {TCs.max():.2f} K\nmin   : {TCs.min():.2f} K'
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)

    # plotting
    plt.figure(figsize=SIZE, dpi=DPI)
    plt.scatter(range(1, len(TCs) + 1), TCs, c='red')
    plt.ylabel('$T_C, K$', fontsize=16)
    plt.grid(alpha=.4)
    plt.text(1, max(TCs), text, verticalalignment='top', bbox=props)
    plt.xticks([])
    plt.savefig(os.path.join(input_folder, 'Tc_plot.pdf'), bbox_inches='tight')


def plot_E_tot(energies_afm: list, nn_matrix: list, input_folder: str):
    x = range(1, len(energies_afm) + 1)
    text = f"""
    $dE$    :  {(energies_afm.max() - energies_afm.min()) * 1000 :.2f}   meV
    max : {energies_afm.max() * 1000:.2f}  meV
    min  : {energies_afm.min() * 1000:.2f}   meV
    """
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)

    # plotting
    plt.figure(figsize=SIZE, dpi=DPI)
    plt.scatter(x, energies_afm)
    plt.grid(alpha=.4)
    plt.xlabel('afm structure id')
    plt.ylabel('$E_{tot}, eV$')
    plt.xticks(x, nn_matrix, rotation=45, ha='right')
    plt.text(1, min(energies_afm), text, verticalalignment='bottom', bbox=props)
    plt.savefig(os.path.join(input_folder, 'Etot_plot.pdf'), bbox_inches='tight')


def plot_E_geom(Eg_list: list, combination_list: list, input_folder: str):
    E_geom = np.array(Eg_list) * 1000
    E_geom_norm = E_geom - E_geom.min()
    max_E_geom = max(E_geom)
    min_E_geom = min(E_geom)
    dE_geom = max_E_geom - min_E_geom
    x_label = range(1, len(E_geom) + 1)
    text = f"""
    $dE$    :  {(max_E_geom - min_E_geom):.2f}   meV
    max : {max_E_geom:.2f}  meV
    min  : {min_E_geom:.2f}   meV
    """
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)

    # plotting
    plt.figure(figsize=SIZE, dpi=DPI)
    plt.scatter(x_label, E_geom_norm)
    plt.title('Normalized energy without exchange coupling')
    plt.ylabel('Eg, meV')
    plt.text(1, max(E_geom_norm), text, verticalalignment='top', bbox=props)
    plt.xticks(x_label, combination_list, rotation=45, ha='right')
    plt.grid(alpha=.4)
    plt.savefig(os.path.join(input_folder, 'Egeom_plot.pdf'), bbox_inches='tight')
