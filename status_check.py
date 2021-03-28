import os
from pymatgen.io.vasp.outputs import Vasprun
import numpy as np
import warnings
from time import sleep
warnings.filterwarnings("ignore")


def vasprun_checker(input_path):
    vasp_inputs_path = os.path.join(input_path, 'vasp_inputs')
    vasprun_pathes = sorted([os.path.join(vasp_inputs_path, i, 'vasprun.xml')
                             for i in os.listdir(vasp_inputs_path)])
    tmp_vasprun = vasprun_pathes.copy()
    while 1:
        print(len(vasprun_pathes))
        for i, vasprun_path in enumerate(vasprun_pathes):
            print(i + 1, end=' ')
            try:
                vasprun = Vasprun(vasprun_path, parse_dos=False,
                                  parse_eigen=False, exception_on_bad_xml=False)
                if vasprun.converged:
                    print(f'Converged! {vasprun_path}')
                    tmp_vasprun.remove(vasprun_path)
                else:
                    print(f'Not converged! {vasprun_path}')
                    tmp_vasprun.remove(vasprun_path)
            except Exception:
                print('Still running')

        vasprun_pathes = tmp_vasprun.copy()
        print('\n')
        sleep(5)
        if not vasprun_pathes:
            print('All done!')
            break


if __name__ == '__main__':
    input_path = input('Enter input path: ')
    vasprun_checker(input_path)
