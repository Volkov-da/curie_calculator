import os
import curie_calculator

number_of_structures = 15
input_folder = 'examples/Bi-Mn/'
curie_calculator.run_enum(input_folder)
curie_calculator.get_structures(num_of_structures=number_of_structures)
curie_calculator.vasp_inputs_creator(num_of_structures=number_of_structures)
curie_calculator.siman_inputs_creator(num_of_structures=number_of_structures)
curie_calculator.enum_out_collector()
