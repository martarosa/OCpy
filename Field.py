import numpy as np
import auxiliary_functions as af

#-------------------------------------------
#Field (x y z) components at time t: 3D vector
#-------------------------------------------



class PropagatorFieldOC():
    def __init__(self):
        self.field_dt = np.array(3)

    def propagate_field_OC_Rabitz(self, ket, bra, muT, alpha):
        f = af.double_summation(ket, np.conj(bra), muT)
        self.field_dt = -1 * np.imag(f) / alpha
        return self.field_dt

    def propagate_field_OC_projector(self, coeff_ket, coeff_bra, muT, alpha):
        f_old = np.dot(np.conj(coeff_bra), np.sum(muT*coeff_ket[np.newaxis, :, np.newaxis], axis=1))
        f = np.dot(np.conj(coeff_ket), coeff_bra) * f_old
        self.field_dt = -1 * np.imag(f) / alpha
        return self.field_dt





#---------------------------------------------
#field matrix at different times to be read or created
#---------------------------------------------
class Field():
    def __init__(self):
        self.field_type = None
        self.dt = None
        self.nstep = None
        self.parameters = dict({})
        self.field = None



    def init_field(self, field_parameters):
        self.field_type = field_parameters.field_type
        self.dt = field_parameters.dt
        self.nstep = field_parameters.nstep
        self.field = field_parameters.field
        self.parameters['fi'] = field_parameters.fi
        self.parameters['fi_cos'] = field_parameters.fi_cos
        self.parameters['omega'] = field_parameters.omega
        self.parameters['sigma'] = field_parameters.sigma
        self.parameters['t0'] = field_parameters.t0
        self.parameters['omega_max'] = field_parameters.omega_max
        self.chose_field(self.field_type)




    def chose_field(self, key):
         field = {
            'const': lambda : self.const_pulse(),
            'pip': lambda : self.pip_pulse(),
            'sin': lambda: self.sin_pulse(),
            'gau': lambda: self.gau_pulse(),
            'sum': lambda: self.sum_pulse(),
            'genetic' : lambda: self.genetic_pulse(),
            'optimizedRabitz' : lambda: self.read_pulse()
             }
         return field.get(key, lambda: "Inexistent field type")()



    def read_pulse(self):
        pass


    def const_pulse(self):
        self.field = np.zeros([self.nstep, 3])
        self.field[:, 0] = self.parameters['fi'][0]
        self.field[:, 1] = self.parameters['fi'][1]
        self.field[:, 2] = self.parameters['fi'][2]




    def pip_pulse(self):
        self.field = np.zeros([self.nstep, 3])
        # Pi pulse: fmax*(cos^2(pi(t-t0)/(2s)) * cos(w(t-t0)))
        for i in range(self.nstep):
            if abs(self.dt * i - self.parameters['t0']) < self.parameters['sigma']:
                self.field[i] = self.parameters['fi'] \
                                * np.square(
                    np.cos(np.pi * (self.dt * i - self.parameters['t0']) / (2 * self.parameters['sigma']))) \
                                * np.cos(self.parameters['omega'] * (self.dt * i - self.parameters['t0']))


    def sin_pulse(self):
        self.field = np.zeros([self.nstep, 3])
        # fmax*exp[-(t-t0)^2/s^2]*sin(wt)
        for i in range(self.nstep):
            self.field[i] = self.parameters['fi'] \
                            * np.exp(-(self.dt * i - self.parameters['t0']) ** 2 / (2 * (self.parameters['sigma']) ** 2)) \
                            * np.sin(self.parameters['omega'] * (self.dt * i))

    def gau_pulse(self):
        self.field = np.zeros([self.nstep, 3])
        for i in range(self.nstep):
            self.field[i] = self.parameters['fi'] \
                            * np.exp(-(self.dt * i - self.parameters['t0']) ** 2 / (self.parameters['sigma']) ** 2)

    def sum_pulse(self):
        self.field = np.zeros([self.nstep, 3])
        # fi*sin(wi*t + fi_cos(wi*t)
        if self.parameters['omega'].ndim == 1:
            for i in range(self.nstep):
                self.field[i] = self.parameters['fi'] * np.sin(self.parameters['omega'] * i * self.dt) \
                                + self.parameters['fi_cos'] * np.cos(self.parameters['omega'] * i * self.dt)
        elif self.parameters['omega'].ndim == 2:
            for i in range(self.nstep):
                for j in range(self.parameters['omega'].shape[0]):
                    self.field[i] = self.field[i] \
                                    + self.parameters['fi'][j] * np.sin(self.parameters['omega'][j] * i * self.dt) \
                                    + self.parameters['fi_cos'][j] * np.cos(self.parameters['omega'][j] * i * self.dt)


    def genetic_pulse(self):
        n_fourier_omega = int(1 + self.parameters['omega_max']*self.dt*self.nstep/(2*np.pi))
        self.parameters['omega'] = np.zeros([n_fourier_omega, 3])
        self.parameters['fi'] = np.full([n_fourier_omega, 3], self.parameters['fi'][0])
        self.parameters['fi_cos'] = np.full([n_fourier_omega, 3], self.parameters['fi'][0])
        for i in range(n_fourier_omega):
            omega_i = np.pi*2*i/(self.dt*self.nstep)
            self.parameters['omega'][i] = [omega_i, omega_i, omega_i]
        self.parameters['sigma'] = 0
        self.parameters['t0'] = 0
        self.sum_pulse()






