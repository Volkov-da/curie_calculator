import os
from curie_calculator import *

number_of_structures = 20
input_folder = 'examples/EuO/'

run_enum(input_folder)
get_structures(num_of_structures=number_of_structures)
vasp_inputs_creator(num_of_structures=number_of_structures)
siman_inputs_creator(num_of_structures=number_of_structures)
enum_out_collector()
submit_jobs(input_folder)
