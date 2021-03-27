from pymatgen.core import Structure
from pymatgen.analysis.magnetism.analyzer import MagneticStructureEnumerator, CollinearMagneticStructureAnalyzer

fm_structure = Structure.from_file('Fe_mp-13_primitive.cif')

fm_structure

mag_prec = 0.1
enum_prec = 0.001

struct_list = MagneticStructureEnumerator(fm_structure, transformation_kwargs={'symm_prec':mag_prec,'enum_precision_parameter':enum_prec})

print(struct_list)