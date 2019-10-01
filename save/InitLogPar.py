import numpy as np

from InitPar import InitPar
from save.LogHeaderParameters import LogHeaderParameters


class InitLogPar(InitPar):
    def __init__(self):
        super().__init__()
        self.parameters = LogHeaderParameters()

    def init(self, user_input):
        self.parameters.dt = user_input.sys.par['dt']
        self.parameters.env = user_input.env.par['env']

        self.parameters.restart = user_input.oc.par['restart']
        self.parameters.target_state = user_input.oc.par["target_state"]
        self.parameters.alpha = user_input.oc.par['alpha']
        self.parameters.alpha0 = user_input.oc.par['alpha0']

        self.parameters.field_type = user_input.field.par['field_type']
        self.parameters.fi = np.array2string(user_input.field.par['fi'])
        self.parameters.fi_cos = np.array2string(user_input.field.par['fi_cos'])
        self.parameters.sigma = user_input.field.par['sigma']
        self.parameters.omega = np.array2string(user_input.field.par['omega'])
        self.parameters.t0 = user_input.field.par['t0']