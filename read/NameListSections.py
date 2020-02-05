import numpy as np
from read.ABCNameListSection import ABCNamelistSection


class SectionSystem(ABCNamelistSection):
    def __init__(self):
        super().__init__()
        self.section = 'SYSTEM'
        self.section_default_dictionary = {'folder': '',
                                           'name': 'output',
                                           'nstep': '10000',
                                           'dt': '0.01',
                                           'oc_algorithm': 'rabitzi'}
        self.section_dictionary = {}
        self.allowed_val = [['oc_algorithm', ['eulero_1order_prop', 'eulero_2order_prop', 'rabitzi', 'rabitzii', 'genetic']]]
        self.case_unsensitive_keys = ['oc_algorithm']
        self.not_implemented_val = []

    def init_default_folder(self, folder):
        self.section_default_dictionary['folder'] = folder



class SectionField(ABCNamelistSection):
    def __init__(self):
        super().__init__()
        self.section = 'FIELD'
        self.section_default_dictionary = {'field_type': 'const',
                                           'fi': '0.01 0.01 0.01',
                                           'omega': '0 0 0',
                                           'sigma': '0',
                                           't0': '0',
                                           'name_field_file': 'false'}
        self.section_dictionary = {}
        self.allowed_val = [['field_type', ['const', 'pip', 'sin', 'gau', 'sum', 'sum_pip', 'genetic', 'test']]]
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
        self.section_default_dictionary = {'name_ci': 'ci_ini.inp',
                                           'name_ei': 'ci_energy.inp',
                                           'name_mut': 'ci_mut.inp'}
        self.section_dictionary = {}
        self.allowed_val = []
        self.case_unsensitive_keys = []
        self.not_implemented_val = []


class SectionEnviron(ABCNamelistSection):
    def __init__(self):
        super().__init__()
        self.section = 'ENVIRON'
        self.section_default_dictionary ={ 'env': 'vac',
                                           'name_vij': 'ci_pot.inp',
                                           'name_q_tdplas': 'np_bem.mdy',
                                           'read_qijn': 'false',
                                           'name_file_qijn': 'qijn.dat',
                                           'name_file_cavity': 'cavity.inp',
                                           'name_q_local_field': 'np_bem.mld'}

        self.section_dictionary = {}
        self.allowed_val = [['env', ['vac', 'sol', 'nanop']]]
        self.case_unsensitive_keys = ['env', 'read_qijn']
        self.not_implemented_val = []


class SectionSave(ABCNamelistSection):
    def __init__(self):
        super().__init__()
        self.section = "SAVE"
        self.section_default_dictionary = {'restart_step': '10'}
        self.section_dictionary = {}
        self.allowed_val = []
        self.case_unsensitive_keys = []
        self.not_implemented_val = []


class SectionOptimalControl(ABCNamelistSection):
    def __init__(self):
        super().__init__()
        self.section = "OPTIMALC"
        self.section_default_dictionary = { 'restart': 'false',
                                            'alpha': 'const',
                                            'alpha0': '1',
                                            'target_state': '1',
                                            'n_iterations': '0',
                                            'convergence_thr': '99999',
                                            'delta_ts': '0',
                                            'Ns': '0',
                                            'iterator_config_file': 'None'}
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



