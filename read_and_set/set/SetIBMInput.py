from read_and_set.set.ABCSetInput import ABCSetInput
from read_and_set.input import IBMInput


class SetIBMInput(ABCSetInput):
    def __init__(self):
        super().__init__()
        self.input_parameters = IBMInput()

    def set(self, user_input):
        self.input_parameters.provider = user_input.parameters.section_dictionary['provider']
        self.input_parameters.device = user_input.parameters.section_dictionary['device']
        self.input_parameters.noise = eval(user_input.parameters.section_dictionary['noise'])
        self.input_parameters.shots = int(user_input.parameters.section_dictionary['shots'])
