import numpy as np
import NameListSections as sec
import configparser

class ReadFieldRestart():
    def __init__(self):
        self.field = sec.SectionField()

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
        self.field = sec.SectionField()


    def read_file(self, folder, namefile):
        user_input = configparser.ConfigParser()
        user_input.read_dict({'field_type': 'optimizedRabitz', #serve per la creazione del campo in field, non fa niente
                              'fi': '0 0 0',
                              'fi_cos':'0 0 0',
                              'omega': '0 0 0',
                              'sigma': '0',
                              't0': '0',
                              'name_field_file': 'false'})

        self.field.par = dict(user_input)
        field = np.loadtxt(folder + namefile, usecols=(0, 1, 2))
        return field