VAMPIRE_PATH = '../../../vampire/vampire-serial'

LDAUJ_DICT = {'Co': 0, 'Cr': 0, 'Fe': 0, 'Mn': 0, 'Mo': 0, 'Ni': 0, 'V': 0, 'W': 0,
              'Nb': 0, 'Sc': 0, 'Ru': 0, 'Rh': 0, 'Pd': 0, 'Cu': 0, 'Y': 0, 'Os': 0, 'Ti': 0, 'Zr': 0, 'Re': 0, 'Hf': 0,
              'Pt': 0, 'La': 0}

LDAUU_DICT = {
    'Co': 3.32, 'Cr': 3.7, 'Fe': 5.3, 'Mn': 3.9, 'Mo': 4.38,
    'Ni': 3.33, 'V': 3.25, 'W': 6.2,
              'Nb': 1.45, 'Sc': 4.18, 'Ru': 4.29, 'Rh': 4.17, 'Pd': 2.96, 'Cu': 7.71, 'Y': 3.23, 'Os': 2.47, 'Ti': 5.89,
              'Zr': 5.55, 'Re': 1.28, 'Hf': 4.77, 'Pt': 2.95, 'La': 5.3}

LDAUL_DICT = {'Co': 2, 'Cr': 2, 'Fe': 2, 'Mn': 2, 'Mo': 2, 'Ni': 2, 'V': 2, 'W': 2,
              'Nb': 2, 'Sc': 2, 'Ru': 2, 'Rh': 2, 'Pd': 2, 'Cu': 2, 'Y': 2, 'Os': 2, 'Ti': 2, 'Zr': 2, 'Re': 2, 'Hf': 2,
              'Pt': 2, 'La': 2}

STAT_DICT = {
    'ENCUT': 550, 'ISMEAR': -5, 'EDIFF': 1E-7, 'SYMPREC': 1E-8, 'NCORE': 4, 'ICHARG': 2,
    'LDAU': True, 'LDAUJ': LDAUJ_DICT, 'LDAUL': LDAUL_DICT, 'LDAUU': LDAUU_DICT, 'NELM': 120,
    'LVHAR': False,
    'LDAUPRINT': 1, 'LDAUTYPE': 2, 'LASPH': True, 'LMAXMIX': 4, 'LWAVE': False, 'LVTOT': False,
    'LAECHG': False
}

DEFAULT_DICT = {
    'ECUT_OPT' : False,
    'KPOINTS_OPT' : False,
    
    'MAGNETIC_ATOM': 'Fe',

    # cut-off energy optimization
    'ECUT_MIN': 400,
    'ECUT_MAX' : 700,
    'ECUT_STEP' : 20,
    
    # Kpoints grid optimization
    'KPOINT_MIN' : 20,
    'KPOINT_MAX' : 80,
    'KPOINT_STEP' : 10,
    
    # monte-carlo simulation vampire
    'MAX_T': 1200,
    'STEP_T' : 20,
    'CUTOFF_RADIUS' : 3.5,
    'NON_MAGNETIC_ATOMS' : 'O',
    'TYPE_OF_CALC' : 'Curie',
    'COUPLING_CONSTANTS' : [1.500e-21, 1.500e-21, 3.381e-21]
    }