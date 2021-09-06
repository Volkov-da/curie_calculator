from pymatgen.core import Structure
from pymatgen.io.vasp.sets import MPStaticSet
from file_builder import create_job_script
import os


def submit_all_jobs(input_path: str, submit_path: str) -> None:
    submit_full_path = os.path.join(input_path, submit_path)
    initial_path = os.getcwd()
    for folder_name in os.listdir(submit_full_path):
        os.chdir(initial_path)
        tmp_path = os.path.join(submit_full_path, folder_name)
        os.chdir(tmp_path)
        os.system('sbatch jobscript.sh')
    os.chdir(initial_path)


def get_ecut_files(input_path: str) -> None:
    ecut_opt_path = os.path.join(input_path, 'encut')
    if not os.path.exists(ecut_opt_path):
        os.mkdir(ecut_opt_path)
    structure = Structure.from_file(os.path.join(input_path, 'POSCAR'))
    for ecut in range(500, 1000, 20):
        user_settings = {'ENCUT': ecut,
                         'EDIFF': 1E-7,
                         'NCORE': 4,
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


if __name__ == '__main__':
    input_path = os.getcwd()
    print(input_path)
    get_ecut_files(input_path)
    submit_all_jobs(input_path=input_path, submit_path='encut')
