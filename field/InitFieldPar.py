import os.path

from field.FieldParameters import FieldParameters
from InitPar import InitPar
from field.ReadFieldRestart import ReadFieldRestart
from read.ReadOutputGaussian import ReadOutputGaussian


class InitFieldPar(InitPar):
    def __init__(self):
        super().__init__()
        self.parameters = FieldParameters()
        self.read_restart = ReadFieldRestart()

    def init(self, user_input):
        if (user_input.oc.par['restart'] == 'false'):
            self.init_no_restart(user_input)
        else:
            self.init_restart(user_input)



    def init_no_restart(self, user_input):
        read_output = ReadOutputGaussian()
        self.parameters.dt = float(user_input.sys.par['dt'])
        self.parameters.nstep = int(user_input.sys.par['nstep'])
        self.parameters.field_type = user_input.field.par['field_type']
        self.parameters.fi = user_input.field.par['fi']
        self.parameters.fi_cos = user_input.field.par['fi_cos']
        self.parameters.sigma = float(user_input.field.par['sigma'])
        self.parameters.omega = user_input.field.par['omega']
        self.parameters.t0 = float(user_input.field.par['t0'])
        self.parameters.omega_max = read_output.read_en_ci0(user_input.sys.par['folder'] +
                                                            user_input.wf.par['name_ei'])[-1]



    def init_restart(self, user_input):
        if os.path.isfile(user_input.sys.par['folder'] + user_input.field.par['name_field_file']):
            if(user_input.sys.par['folder'] + user_input.field.par['name_field_file'] != user_input.sys.par['folder'] + user_input.sys.par['name'] + '_field_bkp.dat'):
                user_input.oc.par['restart'] = 'restart_from_different_name_field'
            self.read_restart.read_file(user_input.sys.par['folder'], user_input.field.par['name_field_file'])
            self.parameters = self.read_restart.field_par
        elif os.path.isfile(user_input.sys.par['folder'] + user_input.sys.par['name'] + '_field_bkp.dat'):
            self.read_restart.read_file(user_input.sys.par['folder'], user_input.sys.par['name'] + '_field_bkp.dat')
            self.parameters = self.read_restart.field_par
        else:
            user_input.oc.par['restart'] = 'norestart_found'
            self.init_no_restart(user_input)