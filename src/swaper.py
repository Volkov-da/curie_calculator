import os
from shutil import copy


def stat_swp(calc_fold: str) -> None:
    log_path = os.path.join(calc_fold, 'log')
    log_path_relax = os.path.join(calc_fold, 'log_relax')
    conctcar_path = os.path.join(calc_fold, 'CONTCAR')
    poscar_path = os.path.join(calc_fold, 'POSCAR')
    incar_relax_path = os.path.join(calc_fold, 'INCAR_relax')
    incar_stat_path = os.path.join(calc_fold, 'INCAR_stat')
    incar_path = os.path.join(calc_fold, 'INCAR')
    copy(log_path, log_path_relax)  # log -> log_relax
    copy(conctcar_path, poscar_path)  # CONTCAR -> POSCAR
    copy(incar_path, incar_relax_path)  # INCAR -> INCAR_relax
    copy(incar_stat_path, incar_path)  # INCAR_stat -> INCAR


def non_conv_swp(calc_fold: str) -> None:
    conctcar_path = os.path.join(calc_fold, 'CONTCAR')
    poscar_path = os.path.join(calc_fold, 'POSCAR')
    copy(conctcar_path, poscar_path)  # CONTCAR -> POSCAR


def swaper(input_path: str) -> None:
    initial_path = os.getcwd()
    converged_msg = ' reached required accuracy - stopping structural energy minimisation'
    swp_msg = 'to POSCAR and continue'
    vasp_inputs_path = os.path.join(input_path, 'vasp_inputs')
    for folder in os.listdir(vasp_inputs_path):
        calc_fold = os.path.join(vasp_inputs_path, folder)
        if os.path.isdir(calc_fold):
            log_path = os.path.join(calc_fold, 'log')
            with open(log_path) as log_f:
                log_text = log_f.readlines()
                last_line = log_text[-1]
                if converged_msg in last_line:
                    stat_swp(calc_fold=calc_fold)
                    print('Converged', calc_fold, end=' ')
                    os.chdir(calc_fold)
                    os.system('sbatch jobscript.sh')
                    os.chdir(initial_path)
                elif swp_msg in last_line:
                    print('Need swap', end=' ')
                    non_conv_swp(calc_fold)
                    os.chdir(calc_fold)
                    os.system('sbatch jobscript.sh')
                    os.chdir(initial_path)


if __name__ == '__main__':
    input_path = os.getcwd()
    swaper(input_path)
