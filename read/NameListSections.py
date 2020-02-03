import numpy as np
from read.ABCNameListSection import ABCNamelistSection



class SectionSystem(ABCNamelistSection):
    def __init__(self):
        super().__init__()
        self.section = 'SYSTEM'
        self.section_keys = ['folder', 'name', 'nstep', 'dt', 'oc_algorithm']
        self.section_default_dictionary = {}
        self.section_dictionary = {}
        self.allowed_val = [['oc_algorithm', ['eulero_1order_prop', 'eulero_2order_prop', 'rabitzi', 'rabitzii', 'genetic']]]
        self.case_unsensitive_keys = ['oc_algorithm']
        self.not_implemented_val = []


class SectionField(ABCNamelistSection):
    def __init__(self):
        super().__init__()
        self.section = 'FIELD'
        self.section_keys = ['field_type', 'fi', 'omega', 'sigma', 't0', 'name_field_file']
        self.section_default_dictionary = {}
        self.section_dictionary = {}
        self.allowed_val = [['field_type', ['const', 'pip', 'sin', 'gau', 'read_file', 'sum', 'sum_pip', 'genetic', 'test']]]
        self.case_unsensitive_keys = ['field_type']
        self.not_implemented_val = []

    def convert_string_coefficients(self, key_string_coeff):
        self.section_dictionary[key_string_coeff] = np.asarray(self.section_dictionary[key_string_coeff].split()).astype(float)
        rows = int(self.section_dictionary[key_string_coeff].shape[0] / 3)
        self.section_dictionary[key_string_coeff].reshape(rows, 3)


class SectionWaveFunction(ABCNamelistSection):
    def __init__(self):
        super().__init__()
        self.section = 'WAVEFUNCTION'
        self.section_keys = ['name_ci','name_ei', 'name_mut']
        self.section_default_dictionary = {}
        self.section_dictionary = {}
        self.allowed_val = []
        self.case_unsensitive_keys = []
        self.not_implemented_val = []


class SectionEnviron(ABCNamelistSection):
    def __init__(self):
        super().__init__()
        self.section = 'ENVIRON'
        self.section_keys = ['env', 'name_vij', 'name_q_tdplas', 'read_qijn', 'name_file_qijn', 'name_file_cavity', 'name_q_local_field']
        self.section_default_dictionary = {}
        self.section_dictionary = {}
        self.allowed_val = [['env', ['vac', 'sol', 'nanop']]]
        self.case_unsensitive_keys = ['env', 'read_qijn']
        self.not_implemented_val = []


class SectionSave(ABCNamelistSection):
    def __init__(self):
        super().__init__()
        self.section = "SAVE"
        self.section_keys = ['restart_step']
        self.section_default_dictionary = {}
        self.section_dictionary = {}
        self.allowed_val = []
        self.case_unsensitive_keys = []
        self.not_implemented_val = []


class SectionOptimalControl(ABCNamelistSection):
    def __init__(self):
        super().__init__()
        self.section = "OPTIMALC"
        self.section_keys = ['restart', 'alpha', 'alpha0', 'target_state', 'n_iterations', 'convergence_thr',
                             'delta_ts', 'Ns', 'iterator_config_file']
        self.section_default_dictionary = {}
        self.section_dictionary = {}
        self.allowed_val = [['alpha', ['const', 'sin', 'quin']],
                            ['restart', ['true', 'false']]]
        self.case_unsensitive_keys = ['restart', 'alpha']


class SectionGenetic(ABCNamelistSection):
    def __init__(self):
        super().__init__()
        self.section = None
        self.section_default_dictionary = {}
        self.section_dictionary = {}
        self.allowed_val = [['mate', ['DEAP_cxUniform']],
                            ['mutate',['DEAP_mutGaussian']],
                            ['select',['DEAP_selBest']],
                            ['deap',['true','false']]]
        self.case_unsensitive_keys = ['deap']


class SectionGenetic_new(ABCNamelistSection):
    def __init__(self):
        super().__init__()
        self.section = None
        self.section_default_dictionary = {}
        self.section_dictionary = {}
        self.allowed_val = [['genetic_algorithm', ['sequential', 'mixed']]]
        self.case_unsensitive_keys = ['genetic_algorithm']


class SectionMate(ABCNamelistSection):
    def __init__(self):
        super().__init__()
        self.section = None
        self.section_default_dictionary = {}
        self.section_dictionary = {}
        self.allowed_val = [['mate', ['DEAP_cxUniform']]]
        self.case_unsensitive_keys = []

class SectionMutate(ABCNamelistSection):
    def __init__(self):
        super().__init__()
        self.section = None
        self.section_default_dictionary = {}
        self.section_dictionary = {}
        self.allowed_val = [['mutate',['DEAP_mutGaussian']]]
        self.case_unsensitive_keys = []

class SectionSelect(ABCNamelistSection):
    def __init__(self):
        super().__init__()
        self.section = None
        self.section_default_dictionary = {}
        self.section_dictionary = {}
        self.allowed_val = [['select',['DEAP_selBest']]]
        self.case_unsensitive_keys = []



