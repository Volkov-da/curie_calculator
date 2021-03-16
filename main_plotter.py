from curie_calculator import get_results, plot_results
import matplotlib.pyplot as plt

FAKE_MAGNETIC_ATOMS = ['Po', 'Mn']
input_folder = 'examples/Bi-Mn/'
num_of_structures = 12
SPIN = 5/2

out_dict = get_results(input_folder, num_of_structures, FAKE_MAGNETIC_ATOMS, SPIN)

Tcs = out_dict['Tc_list']

plot_results(Tcs)
