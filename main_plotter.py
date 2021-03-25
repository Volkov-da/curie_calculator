from curie_calculator import input_reader
from solver import get_results
from plotter import plot_E_geom, plot_E_tot, plot_TCs

input_folder, num_of_structures, magnetic_atoms, spin = input_reader()

print(f'Working with "{num_of_structures}" structures from {input_folder} folder:')

out_dict = get_results(input_folder, num_of_structures, magnetic_atoms, spin)
plot_E_tot(out_dict['energies_afm'], out_dict['nn_matrix'], input_folder)
plot_TCs(out_dict['Tc_list'], input_folder)
plot_E_geom(out_dict['Eg_list'], out_dict['combination_list'], input_folder)

print(out_dict)

print(f'Success! Now you can check your "{input_folder}" folder for results.')
