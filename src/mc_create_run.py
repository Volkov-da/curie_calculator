import os
import numpy as np
from read_input import update_defaults
from solver import get_spin
from variables import DEFAULT_DICT, VAMPIRE_PATH
from pymatgen.core.structure import Structure


def input_file_vamp(input_path: str, t_max: int, t_step: int) -> None:
    out_path = os.path.join(input_path, 'monte_carlo')
    if not os.path.exists(out_path):
        os.mkdir(out_path)
    input_text = f"""#------------------------------------------
# Creation attributes:
#------------------------------------------
create:crystal-structure=sc
#------------------------------------------
# System Dimensions:
#------------------------------------------
dimensions:system-size-x = 5 !nm
dimensions:system-size-y = 5 !nm
dimensions:system-size-z = 5 !nm
#------------------------------------------
# Material Files:
#------------------------------------------
material:file=structure.mat
material:unit-cell-file=structure.ucf
#------------------------------------------
# Simulation attributes:
#------------------------------------------
sim:minimum-temperature=0.0
sim:maximum-temperature={t_max}
sim:temperature-increment={t_step}

sim:time-steps-increment=100
sim:equilibration-time-steps=10000
sim:loop-time-steps=30000

sim:time-step=1.0E-16
#------------------------------------------
# Program and integrator details
#------------------------------------------
sim:program=curie-temperature
sim:integrator=monte-carlo
#------------------------------------------
# data output
#------------------------------------------
output:temperature
output:mean-magnetisation-length
output:material-mean-magnetisation-length

screen:temperature
screen:mean-magnetisation-length
screen:material-mean-magnetisation-length"""
    file_out_path = os.path.join(out_path, 'input')
    with open(file_out_path, 'w+') as out_f:
        out_f.writelines(input_text)


def mat_file_vamp(input_path: str, magmom: float, magnetic_atom: str, type_of_calc: str):
    num_materials = 1 if type_of_calc == 'Curie' else 2
    input_text = f"""#---------------------------------------------------
# Number of Materials
#---------------------------------------------------
material:num-materials={num_materials}
#---------------------------------------------------
# Material 1 {magnetic_atom} Spin-UP
#---------------------------------------------------
material[1]:material-name={magnetic_atom}-up
material[1]:damping-constant=1.0
material[1]:atomic-spin-moment={magmom} !muB
material[1]:material-element={magnetic_atom}
material[1]:uniaxial-anisotropy-constant=0.0
material[1]:initial-spin-direction=0,0,1"""
    if type_of_calc == 'Neel':
        input_text += f"""
#---------------------------------------------------
# Material 2 {magnetic_atom} Spin-DOWN
#---------------------------------------------------
material[2]:material-name=Fe-dn
material[2]:damping-constant=1.0
material[2]:atomic-spin-moment={magmom} !muB
material[2]:material-element={magnetic_atom}
material[2]:uniaxial-anisotropy-constant=0.0
material[2]:initial-spin-direction=0,0,-1"""
    out_path = os.path.join(input_path, 'monte_carlo')
    if not os.path.exists(out_path):
        os.mkdir(out_path)

    file_out_path = os.path.join(out_path, 'structure.mat')
    with open(file_out_path, 'w+') as out_f:
        out_f.writelines(input_text)


def ucf_file_vamp(input_path: str,  coupling_constants: list, non_magnetic_atoms: list, type_of_calc: str, cutoff_radius):
    # set 1 for Curie temperature, 2 for Néel temperature
    num_materials = 1 if type_of_calc == 'Curie' else 2
    poscar_path = os.path.join(input_path, 'POSCAR')
    structure = Structure.from_file(poscar_path)
    structure.remove_species(non_magnetic_atoms)
    center_indices, point_indices, offset_vectors, distances = structure.get_neighbor_list(
        cutoff_radius)
    DISTANCES_BETWEEN_NEIGHBORS = sorted(np.unique(distances.round(3)))

    # get thresholds for comparing atomic distances
    thresholds = [(a + b) / 2 for a, b in zip(DISTANCES_BETWEEN_NEIGHBORS[:-1],
                                              DISTANCES_BETWEEN_NEIGHBORS[1:])]
    thresholds = [0.0] + thresholds + [100.0]

    # print unit cell size for VAMPIRE input
    out_path = os.path.join(input_path, 'monte_carlo')
    if not os.path.exists(out_path):
        os.mkdir(out_path)
    out_file_path = os.path.join(out_path, 'structure.ucf')
    with open(out_file_path, 'w') as f:
        f.write('# Unit cell size:\n')
        f.write(
            f'{structure.lattice.a:.6f}  {structure.lattice.b:.6f}  {structure.lattice.c:.6f}\n')

        f.write('# Unit cell vectors:\n')
        for i, (vector, modulus) in enumerate(zip(structure.lattice.matrix, structure.lattice.abc)):
            v = vector / modulus
            f.write(f'{v[0]: .6f}   {v[1]: .6f}   {v[2]: .6f}\n')

        # print fractional coordinates for VAMPIRE input
        f.write('# Atoms num_atoms num_materials; id cx cy cz mat cat hcat\n')
        f.write(f'{len(structure):d} {num_materials}\n')
        for i, coord in enumerate(structure.frac_coords):
            for j, item in enumerate(coord):
                if item < 0:
                    coord[j] += 1
                elif item >= 1:
                    coord[j] -= 1
            material = i // (len(structure) // num_materials)
            f.write(
                f'{i:2d}   {coord[0]:.6f}  {coord[1]:.6f}  {coord[2]:.6f}  {material} 0 0\n')

        f.write('# Interactions n exctype; id i j dx dy dz Jij\n')
        f.write(f'{len(center_indices)} isotropic\n')
        for i, (atom1, atom2, offset, distance) in enumerate(zip(center_indices, point_indices, offset_vectors, distances)):
            for low_lim, high_lim, coupling_constant in zip(thresholds[:-1], thresholds[1:], coupling_constants):
                if low_lim < distance < high_lim:
                    f.write(f'{i:3d}   {atom1:2d}  {atom2:2d}  '
                            f'{int(offset[0]):2d} {int(offset[1]):2d} {int(offset[2]):2d}   {coupling_constant: 6.4e}\n')


def job_monte_carlo(input_path: str, vampire_path: str, job_id=None) -> None:
    if not job_id:
        job_id = os.path.basename(input_path)

    out_path = os.path.join(input_path, 'monte_carlo')
    job_script_text = f"""#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=12:00:00
#SBATCH --job-name={job_id}_monte_carlo
#SBATCH --output=log
#SBATCH --error=err
#SBATCH -p lenovo
{vampire_path}"""
    if not os.path.exists(input_path):
        os.mkdir(input_path)
    file_out_path = os.path.join(out_path, 'jobscript.sh')
    with open(file_out_path, 'w') as job:
        job.writelines(job_script_text)


def submit_monte_carlo(input_path: str) -> None:
    initial_path = os.getcwd()
    tmp_path = os.path.join(input_path, 'monte_carlo')
    os.chdir(tmp_path)
    os.system('sbatch jobscript.sh')
    os.chdir(initial_path)

def run_monte_carlo(input_path: str):
    update_defaults(input_path, DEFAULT_DICT)
    print(DEFAULT_DICT)
    input_file_vamp(
        input_path,
        t_max  = DEFAULT_DICT['MAX_T'],
        t_step = DEFAULT_DICT['STEP_T']

    )

    fm_outcar_path = os.path.join(input_path, 'vasp_inputs', 'fm0', 'OUTCAR')
    ref_magmom = get_spin(fm_outcar_path)

    mat_file_vamp(
        input_path,
        magmom        = ref_magmom,
        magnetic_atom = DEFAULT_DICT['MAGNETIC_ATOM'],
        type_of_calc  = DEFAULT_DICT['TYPE_OF_CALC']
    )

    ucf_file_vamp(
        input_path,
        coupling_constants = DEFAULT_DICT['COUPLING_CONSTANTS'],
        non_magnetic_atoms = DEFAULT_DICT['NON_MAGNETIC_ATOMS'],
        type_of_calc       = DEFAULT_DICT['TYPE_OF_CALC'],
        cutoff_radius      = DEFAULT_DICT['CUTOFF_RADIUS']
    )

    job_monte_carlo(
        input_path,
        vampire_path=VAMPIRE_PATH)

    submit_monte_carlo(input_path)

if __name__ == '__main__':
    run_monte_carlo(input_path=os.getcwd())
