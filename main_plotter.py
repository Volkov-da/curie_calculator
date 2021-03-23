from curie_calculator import get_results, plot_TCs, plot_Etot, input_reader
import matplotlib.pyplot as plt

input_folder, num_of_structures, fake_magnetic_atoms, spin = input_reader()
out_dict = get_results(input_folder, num_of_structures, fake_magnetic_atoms, spin)
plot_TCs(out_dict['Tc_list'], input_folder=input_folder)
plot_Etot(out_dict['energies_afm'], out_dict['nn_matrix'], input_folder=input_folder)
print(out_dict)
