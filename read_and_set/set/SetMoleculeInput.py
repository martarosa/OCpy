from read_and_set.input.MoleculeInput import MoleculeInput
from read_and_set.read.external_output.ReadOutputGaussian import ReadOutputGaussian
from read_and_set.set.ABCSetInput import ABCSetInput


class SetMoleculeInput(ABCSetInput):
    def __init__(self):
        super().__init__()
        self.input_parameters = MoleculeInput()


    def set(self, user_input):
        read_output = ReadOutputGaussian()
        self.input_parameters.wf_ci = read_output.read_ci0(user_input.sys.section_dictionary['folder'] +
                                                           user_input.wf.section_dictionary['name_ci'])

        self.input_parameters.en_ci = read_output.read_en_ci0(user_input.sys.section_dictionary['folder'] +
                                                              user_input.wf.section_dictionary['name_ei'])
        self.input_parameters.muT = read_output.read_muT(user_input.sys.section_dictionary['folder'] +
                                                         user_input.wf.section_dictionary['name_mut'], self.input_parameters.wf_ci.size)
        if user_input.medium.section_dictionary['medium'] != 'vac':
            self.input_parameters.Vijn = read_output.read_V(user_input.sys.section_dictionary['folder'] +
                                                            user_input.medium.section_dictionary['name_vij'],
                                                            self.input_parameters.wf_ci.size)