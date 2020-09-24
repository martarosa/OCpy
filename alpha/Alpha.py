import numpy as np

from parameters.AlphaParameters import AlphaParameters
from SystemObj import DiscreteTimePar


class Alpha():
    def __init__(self):
        self.alpha_t = None
        self.par = AlphaParameters()

    def init_alpha(self, alpha_input):
        self.par.alpha_type = alpha_input.alpha_type
        self.par.read_name = alpha_input.read_name


        discrete_t_par = DiscreteTimePar()
        discrete_t_par.dt = alpha_input.dt
        discrete_t_par.nstep = alpha_input.nstep

        alpha = alpha_input.alpha


        self.chose_alpha(self.par.alpha_type, discrete_t_par, alpha)


    def chose_alpha(self, key, discrete_t_par = None, alpha = None):
        alpha_choice = {
            'const': lambda: self.const_alpha(discrete_t_par),
            'sin': lambda: self.sin_alpha(discrete_t_par),
            'quin': lambda: self.quin_alpha(discrete_t_par),
            'read': lambda: self.read_alpha()
        }
        return alpha_choice.get(key, lambda: "Inexistent alpha type")()

    def read_alpha(self):
        self.alpha_t = np.loadtxt(self.par.read_name)


    def const_alpha(self, discrete_t_par):
        self.alpha_t = np.zeros([discrete_t_par.nstep])
        for i in range(discrete_t_par.nstep):
            self.alpha_t[i] = 1

    def sin_alpha(self, discrete_t_par):
        self.alpha_t = np.zeros([discrete_t_par.nstep])
        for i in range(discrete_t_par.nstep):
            self.alpha_t[i] = 1 / np.square(np.sin(np.pi * (i + 1) / discrete_t_par.nstep))

    def quin_alpha(self, discrete_t_par):
        self.alpha_t = np.zeros([discrete_t_par.nstep])
        for i in range(discrete_t_par.nstep):
                self.alpha_t[i] = 1 / np.exp(
                    -np.power((((i) * discrete_t_par.dt - 125) / 220), 12))