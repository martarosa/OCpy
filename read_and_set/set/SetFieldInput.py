from read_and_set.input.FieldInput import FieldInput
from read_and_set.set.ABCSetInput import ABCSetInput
from read_and_set.read.ReadOutputQuantumCalc import ReadOutputQuantumCalc


class SetFieldInput(ABCSetInput):
    def __init__(self):
        super().__init__()
        self.input_parameters = FieldInput()
        self.read_restart = None #ABCReadFieldRestart()

    def set(self, user_input):
        if (user_input.oc.section_dictionary['restart'] == 'false' or user_input.oc.section_dictionary['restart'] == 'norestart_found'):
            self.set_no_restart(user_input)
        else:
            self.set_restart(user_input)

    def set_no_restart(self, user_input):
        read_output = ReadOutputQuantumCalc()
        self.input_parameters.dt = float(user_input.sys.section_dictionary['dt'])
        self.input_parameters.nstep = int(user_input.field.section_dictionary['nstep'])
        self.input_parameters.additional_steps = int(user_input.field.section_dictionary['additional_steps'])
        self.input_parameters.field_type = user_input.field.section_dictionary['field_type']
        self.input_parameters.fi = user_input.field.section_dictionary['fi']
        self.input_parameters.sigma = float(user_input.field.section_dictionary['sigma'])
        self.input_parameters.omega = user_input.field.section_dictionary['omega']
        self.input_parameters.t0 = float(user_input.field.section_dictionary['t0'])

        self.input_parameters.omega_sys = read_output.read_en_ci0(user_input.sys.section_dictionary['folder'] +
                                                                  user_input.wf.section_dictionary['name_ei'])
        self.input_parameters.field_restart_name = user_input.sys.section_dictionary['folder'] + user_input.field.section_dictionary['name_field_file']

    def set_restart(self, user_input):
            self.read_restart.read_file(user_input.sys.section_dictionary['folder'], user_input.field.section_dictionary['name_field_file'])
            if(user_input.sys.section_dictionary['oc_algorithm'] == 'genetic'):
                self.input_parameters = self.read_restart.input_parameters
                self.input_parameters.dt = float(user_input.sys.section_dictionary['dt'])
                self.input_parameters.nstep = int(user_input.field.section_dictionary['nstep'])
                self.input_parameters.additional_steps = int(user_input.field.section_dictionary['additional_steps'])
            else:
                self.input_parameters = self.read_restart.input_parameters

