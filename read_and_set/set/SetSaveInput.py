from read_and_set.set.ABCSetInput import ABCSetInput
from read_and_set.input.SaveInput import SaveInput


class SetSaveInput(ABCSetInput):
    def __init__(self):
        super().__init__()
        self.input_parameters = SaveInput()

    def set(self, user_input):
        self.input_parameters.folder = user_input.sys.section_dictionary['folder']
        self.input_parameters.name = user_input.sys.section_dictionary['name']
        self.input_parameters.restart_calculation = user_input.oc.section_dictionary['restart']
        self.input_parameters.restart_step = int(user_input.save.section_dictionary['restart_step'])
        self.input_parameters.append = user_input.save.section_dictionary['append']