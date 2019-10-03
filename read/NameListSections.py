import sys
import numpy as np

class Section():
    def __init__(self):
        self.par = None
        self.allowed_val = None
        self.case_unsensitive_keys = None

    def lowering_case(self):
        for key in self.case_unsensitive_keys:
            self.par[key] = self.par[key].lower()

    def check_parameters_allowed_values(self):
        if self.allowed_val:
            for key, par in self.allowed_val:
                if self.par[key] not in par:
                    sys.exit("Error in input. Wrong \"" + key + "\" value")


class SectionSystem(Section):
    def __init__(self):
        super().__init__()
        self.allowed_val = [['propagation', ['eulero_1order', 'eulero_2order', 'rabitzi', 'rabitzii', 'genetic']]]
        self.case_unsensitive_keys = ['propagation']


class SectionField(Section):
    def __init__(self):
        super().__init__()
        self.allowed_val = [['field_type', ['const', 'pip', 'sin', 'gau', 'read_file', 'sum', 'optimizedRabitz', 'genetic']]]
        self.case_unsensitive_keys = ['field_type']

    def convert_string_coefficients(self, key_string_coeff):
        self.par[key_string_coeff] = np.asarray(self.par[key_string_coeff].split()).astype(float)
        rows = int(self.par[key_string_coeff].shape[0]/3)
        self.par[key_string_coeff].reshape(rows, 3)



class SectionWaveFunction(Section):
    def __init__(self):
        super().__init__()
        self.allowed_val = []
        self.case_unsensitive_keys = []


class SectionEnviron(Section):
    def __init__(self):
        super().__init__()
        self.allowed_val = [['env', ['vac', 'sol', 'nanop']]]
        self.case_unsensitive_keys = ['env', 'read_qijn']

class SectionSave(Section):
    def __init__(self):
        super().__init__()
        self.allowed_val = []
        self.case_unsensitive_keys = []


class SectionOptimalControl(Section):
    def __init__(self):
        super().__init__()
        self.allowed_val = [['alpha', ['const', 'sin', 'quin']],
                            ['restart', ['true', 'false']]]
        self.case_unsensitive_keys = ['restart', 'alpha']

class SectionGenetic(Section):
    def __init__(self):
        super().__init__()
        self.allowed_val = [['mate', ['cxUniform']],
                            ['mutate',['mutGaussian']],
                            ['select',['selBest']]]
        self.case_unsensitive_keys = []