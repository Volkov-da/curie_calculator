from curie_calculator import get_results, plot_results, input_reader
import matplotlib.pyplot as plt

fake_magnetic_atoms, input_folder, num_of_structures, spin = input_reader()
out_dict = get_results(input_folder, num_of_structures, fake_magnetic_atoms, spin)
Tcs = out_dict['Tc_list']
plot_results(Tcs, input_folder)
print(out_dict)
