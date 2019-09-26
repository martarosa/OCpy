import numpy as np
import NameListSections as sec
from FieldParameters import FieldParameters
import configparser

class ReadFieldRestart():
    def __init__(self):
        self.field_par = sec.SectionField()

    def read_file(self, folder, namefile):
        pass


class ReadFieldRestartGenetic(ReadFieldRestart):
    def __init__(self):
        super().__init__()
        self.field = sec.SectionField()


    def read_file(self, folder, namefile):
        user_input = configparser.ConfigParser()
        user_input.read_dict({'field_type': 'const',
                              'fi': '0 0 0',
                              'fi_cos':'0 0 0',
                              'omega': '0 0 0',
                              'sigma': '0',
                              't0': '0',
                              'name_field_file': 'false'})

        user_input.read(folder + namefile)
        self.field.par = dict(user_input)
        self.field.lowering_case()
        self.field.check_parameters_allowed_values()
        self.field.convert_string_coefficients('fi')
        self.field.convert_string_coefficients('fi_cos')
        self.field.convert_string_coefficients('omega')
        field = None
        return field



class ReadFieldRestartRabitz(ReadFieldRestart):
    def __init__(self):
        super().__init__()
        self.field_par = FieldParameters


    def read_file(self, folder, namefile):
        self.field_par.field = np.loadtxt(folder + namefile, usecols=(0, 1, 2))
        self.field_par. field_type = 'optimizedRabitz'
        self.field_par.fi = np.array([0, 0, 0])
        self.field_par.fi_cos = np.array([0, 0, 0])
        self.field_par.omega = np.array([0, 0, 0])
        self.field_par.sigma = 0
        self.field_par.t0 = 0
        self.field_par.namefile = 'false'
        self.field_par.omega_max = 0
