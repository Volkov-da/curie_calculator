from stat_file_builder import file_builder
from solver import solver
from time import gmtime, strftime
from os import getcwd

if __name__ == '__main__':
    input_path = getcwd()
    magnetic_atom = input('Enter magnetic atom (str): ')
    print(input_path)
    file_builder(input_path)
    print(strftime("%H:%M:%S", gmtime()), 'Energy calculations done! Starting Hamiltonian fitting\n')
    solver(input_path, magnetic_atom)
