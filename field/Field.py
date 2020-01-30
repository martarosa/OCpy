import numpy as np
from read import auxiliary_functions as af
import random

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



    def init_field(self, field_input):
        self.field_type = field_input.field_type
        self.dt = field_input.dt
        self.nstep = field_input.nstep
        self.field = field_input.field
        self.parameters['fi'] = field_input.fi
        self.parameters['omega'] = field_input.omega
        self.parameters['sigma'] = field_input.sigma
        self.parameters['t0'] = field_input.t0
        self.parameters['omega_sys'] = field_input.omega_sys
        self.chose_field(self.field_type)




    def chose_field(self, key):
         field = {
            'const': lambda : self.const_pulse(),
            'pip': lambda : self.pip_pulse(),
            'sin': lambda: self.sin_pulse(),
            'gau': lambda: self.gau_pulse(),
            'sum': lambda: self.sum_pulse(),
            'sum_pip': lambda: self.sum_pip_pulse(),
            'genetic' : lambda: self.test_pulse(),
            'test'    :lambda: self.test_fourier_pulse(),
             #only internal values
            'restart_rabitz' : lambda: self.restart_rabitz(),
            'restart_genetic' : lambda: self.sum_pulse()
             }
         return field.get(key, lambda: "Inexistent field type")()



    def restart_rabitz(self):
        #everithing is already done in restart
        pass


    def const_pulse(self):
        self.field = np.zeros([self.nstep, 3])
        self.field[:, 0] = self.parameters['fi'][0]
        self.field[:, 1] = self.parameters['fi'][1]
        self.field[:, 2] = self.parameters['fi'][2]


    def test_fourier_pulse(self):
        self.parameters['omega'] = np.array(
            [[0, 0, 0], [0.1266749, 0.1266749, 0.1266749], [0.32140573, 0.32140573, 0.32140573]])
        a=np.array([[0.1,0,0],[0,0.02,0],[0,0,0.1]])
        self.parameters['fi'] = a
        self.sum_pulse()


    def test_pulse(self):
        a = []
        for i in range(15):
            a.append(round(random.uniform(-0.005, 0.005),4))
        a=np.array(a)
        a=a.reshape((5,3))
        self.parameters['omega'] = np.array([[0,0,0],[0.1266749, 0.1266749, 0.1266749],[0.32140573, 0.32140573, 0.32140573], [0.33407322, 0.33407322,0.33407322],[0.34652756,0.34652756,0.34652756]])
        self.parameters['sigma'] = 125
        self.parameters['t0'] = 125
        self.parameters['fi'] = a
        self.sum_pip_pulse()


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
        # f0+fi*sin(wi*t)
        if self.parameters['omega'].ndim == 1: #only one value of omega, 3D field [wx,wy,wz]
            for i in range(self.nstep):
                self.field[i] = self.parameters['fi'][0] - self.parameters['fi'][0] * np.sin(self.parameters['omega'] * i * self.dt) \
                                + self.parameters['fi'] * np.sin(self.parameters['omega'] * i * self.dt)
        elif self.parameters['omega'].ndim == 2: #matrix, N values of omega, 3D [w0x,w0y,w0z],[w1x,w1y,w1z]... il primo valore di omega deve essere sempre 0
            for i in range(self.nstep):
                self.field[i] = self.parameters['fi'][0]
                for j in range(self.parameters['omega'].shape[0]):
                    self.field[i] = self.field[i] \
                                    + self.parameters['fi'][j] * np.sin(self.parameters['omega'][j] * i * self.dt)



    def sum_pip_pulse(self):
        self.field = np.zeros([self.nstep, 3])
        # f0+fi*sin(wi*t)
        for i in range(self.nstep):
            self.field[i] = self.parameters['fi'][0]
            if abs(self.dt * i - self.parameters['t0']) < self.parameters['sigma']:
                for j in range(self.parameters['omega'].shape[0]):
                    self.field[i] = self.field[i] \
                                    + self.parameters['fi'][j] \
                                    * np.square( np.cos(np.pi * (self.dt * i - self.parameters['t0'])/(2 * self.parameters['sigma']))) \
                                    * np.cos(self.parameters['omega'][j] * (i * self.dt - self.parameters['t0']))

    def genetic_pulse(self):
        self.chose_omega_fourier()
        print(self.parameters['omega'].shape)
        self.parameters['fi'] = np.full([self.parameters['omega'].shape[0]+1, 3], self.parameters['fi'][0])
        self.sum_pulse()


    def chose_omega_fourier(self):
        n_fourier_omega_min = 0
        n_fourier_omega_max = int(2 + self.parameters['omega_sys'][-1]*self.dt*self.nstep/(2*np.pi))
        n_fourier_omega = n_fourier_omega_max - n_fourier_omega_min
        self.parameters['omega'] = np.zeros([n_fourier_omega, 3])
        for n in range(n_fourier_omega):
            omega_n = (np.pi*2*(n + n_fourier_omega_min))/(self.dt*self.nstep)
            self.parameters['omega'][n] = [omega_n, omega_n, omega_n]


    def chose_omega_energy(self):
        freq= []
        for i in range(self.parameters['omega_sys'].shape[0]):
            for j in range(i):
                freq.append(self.parameters['omega_sys'][i]-self.parameters['omega_sys'][j])
        self.parameters['omega'] = np.zeros([len(freq), 3])
        self.parameters['omega'][:,0] = freq
        self.parameters['omega'][:, 1] = freq
        self.parameters['omega'][:,2] = freq
