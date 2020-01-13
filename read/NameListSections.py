import sys
import numpy as np



class Section():
    def __init__(self):
        self.section = None
        self.section_default_dictionary = {}
        self.section_dictionary = {}
        self.allowed_val = None
        self.case_unsensitive_keys = None


    def init(self, user_input, default, section):
        self.section_dictionary = dict(user_input[section])
        self.section_default_dictionary = default[section]


    def check(self):
        self.check_keys()
        self.lowering_case()
        self.check_allowed_values()



    def check_keys(self):
            if not all(elem in self.section_dictionary.keys() for elem in self.section_default_dictionary.keys()):
                sys.exit("Error in input. Wrong key in " + self.section + "\n")


    def lowering_case(self):
        for key in self.case_unsensitive_keys:
            self.section_dictionary[key] = self.section_dictionary[key].lower()


    def check_allowed_values(self):
        if self.allowed_val:
            for key, par in self.allowed_val:
                if self.section_dictionary[key] not in par:
                    sys.exit("Error in input. Wrong \"" + key + "\" value. \n"
                             + "Default value is: " + self.section_default_dictionary[key] + "\n")







class SectionSystem(Section):
    def __init__(self):
        super().__init__()
        self.section = 'SYSTEM'
        self.section_default_dictionary = {}
        self.section_dictionary = {}
        self.allowed_val = [['propagation', ['eulero_1order', 'eulero_2order', 'rabitzi', 'rabitzii', 'genetic']]]
        self.case_unsensitive_keys = ['propagation']
        self.not_implemented_val = []


class SectionField(Section):
    def __init__(self):
        super().__init__()
        self.section = 'FIELD'
        self.section_default_dictionary = {}
        self.section_dictionary = {}
        self.allowed_val = [['field_type', ['const', 'pip', 'sin', 'gau', 'read_file', 'sum', 'optimizedRabitz', 'genetic']]]
        self.case_unsensitive_keys = ['field_type']
        self.not_implemented_val = []

    def convert_string_coefficients(self, key_string_coeff):
        self.section_dictionary[key_string_coeff] = np.asarray(self.section_dictionary[key_string_coeff].split()).astype(float)
        rows = int(self.section_dictionary[key_string_coeff].shape[0] / 3)
        self.section_dictionary[key_string_coeff].reshape(rows, 3)


class SectionWaveFunction(Section):
    def __init__(self):
        super().__init__()
        self.section = 'WAVEFUNCTION'
        self.section_default_dictionary = {}
        self.section_dictionary = {}
        self.allowed_val = []
        self.case_unsensitive_keys = []
        self.not_implemented_val = []


class SectionEnviron(Section):
    def __init__(self):
        super().__init__()
        self.section = 'ENVIRON'
        self.section_default_dictionary = {}
        self.section_dictionary = {}
        self.allowed_val = [['env', ['vac', 'sol', 'nanop']]]
        self.case_unsensitive_keys = ['env', 'read_qijn']
        self.not_implemented_val = []


class SectionSave(Section):
    def __init__(self):
        super().__init__()
        self.section = "SAVE"
        self.section_default_dictionary = {}
        self.section_dictionary = {}
        self.allowed_val = []
        self.case_unsensitive_keys = []
        self.not_implemented_val = []

class SectionOptimalControl(Section):
    def __init__(self):
        super().__init__()
        self.section = "OPTIMALC"
        self.section_default_dictionary = {}
        self.section_dictionary = {}
        self.allowed_val = [['alpha', ['const', 'sin', 'quin']],
                            ['restart', ['true', 'false']]]
        self.case_unsensitive_keys = ['restart', 'alpha']


class SectionGenetic(Section):
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


class SectionTDPlasPropagate(Section):
    def __init__(self):
        super().__init__()
        self.section = 'propagate'
        self.section_default_dictionary = {}
        self.section_dictionary = {}
        self.allowed_val = [['out_level', ['high', 'low']],
                            ['interaction_type',['pcm']],
                            ['propagation_type',['ief']],
                            ['interaction_init', ['non-scf']],
                            ['debug_type', ['non']],
                            ['medium_relax', ['non']],
                            ['test_type',['non']]]
        self.case_unsensitive_keys = ['out_level',
                                      'interaction_type',
                                      'propagation_type',
                                      'interaction_init',
                                      'debug_type',
                                      'medium_relax',
                                      'test_type']



class SectionTDPlasMedium(Section):
    def __init__(self):
        super().__init__()
        self.section = 'medium'
        self.section_default_dictionary = {}
        self.section_dictionary = {}
        self.allowed_val = [['medium_type', ['nan', 'sol']],
                            ['medium_init',['fro']],
                            ['medium_pol',['chr']],
                            ['bem_type', ['diag']],
                            ['bem_read_write', ['rea']]]
        self.case_unsensitive_keys = ['medium_type',
                                      'medium_init',
                                      'medium_pol',
                                      'bem_type',
                                      'bem_read_write']


class SectionTDPlasSurface(Section):
    def __init__(self):
        super().__init__()


class SectionTDPlasEps(Section):
    def __init__(self):
        super().__init__()
        self.section = 'eps_function'
        self.section_default_dictionary = {}
        self.section_dictionary = {}
        self.allowed_val = [['epsilon_omega', ['drl']]]
        self.case_unsensitive_keys = ['epsilon_omega']
