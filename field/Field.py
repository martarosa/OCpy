import numpy as np

from SystemObj import DiscreteTimePar, Func_tMatrix

#---------------------------------------------
#field matrix at different times to be read or created
#---------------------------------------------
from parameters.FieldParameters import FieldParameters


class Field():
    def __init__(self):
        self.field = Func_tMatrix()
        self.par = FieldParameters()


    def init_field(self, field_input):
        self.par.field_type = field_input.field_type
        self.par.fi = field_input.fi
        self.par.omega = field_input.omega
        self.par.sigma = field_input.sigma
        self.par.t0 = field_input.t0
        self.par.omega_sys = field_input.omega_sys
        self.par.restart_name = field_input.field_restart_name

        discrete_t_par = DiscreteTimePar()
        discrete_t_par.dt = field_input.dt
        discrete_t_par.nstep = field_input.nstep

        field = Func_tMatrix()
        field.f_xyz = field_input.field
        field.time_axis = field_input.field_time_axis
        self.chose_field(self.par.field_type, discrete_t_par, field)



    def chose_field(self, key, discrete_t_par = None, field = None):
         field_choice = {
            'const': lambda : self.const_pulse(discrete_t_par),
            'pip': lambda : self.pip_pulse(discrete_t_par),
            'sin': lambda: self.sin_pulse(discrete_t_par),
            'gau': lambda: self.gau_pulse(discrete_t_par),
            'sum': lambda: self.sum_pulse(discrete_t_par),
            'genetic' : lambda: self.genetic_pulse(discrete_t_par),
            'sin_cos': lambda : self.sin_cos_pulse(discrete_t_par),
            'read': lambda : self.read_pulse(discrete_t_par),
            'read_genetic': lambda : self.read_genetic_pulse(discrete_t_par),
             #only internal values
            'restart_rabitz' : lambda: self.restart_rabitz(field),
            'restart_genetic' : lambda: self.restart_genetic(discrete_t_par)
             }
         return field_choice.get(key, lambda: "Inexistent field type")()




    def restart_rabitz(self, field):
        self.field = field

    def read_pulse(self, discrete_t_par):
        self.field.time_axis = np.loadtxt(self.par.restart_name, usecols=(2))
        self.field.f_xyz = np.loadtxt(self.par.restart_name, usecols = (3,4,5))

    def read_genetic_pulse(self, discrete_t_par):
        load = np.loadtxt(self.par.restart_name) # w0, w1, ...wn \n a0 a1 ...an
        n = int((load.shape[0])/2)
        self.par.fi = load[n:]
        print(self.par.fi)
        self.par.omega = load[:n]
        self.sum_pulse(discrete_t_par)



    def const_pulse(self, discrete_t_par):
        self.field.time_axis = np.linspace(0, discrete_t_par.dt * (discrete_t_par.nstep-1), discrete_t_par.nstep)
        self.field.f_xyz = np.zeros([discrete_t_par.nstep, 3])
        self.field.f_xyz[:, 0] = self.par.fi[0]
        self.field.f_xyz[:, 1] = self.par.fi[1]
        self.field.f_xyz[:, 2] = self.par.fi[2]




    def pip_pulse(self, discrete_t_par):
        self.field.time_axis = np.linspace(0, discrete_t_par.dt * (discrete_t_par.nstep - 1), discrete_t_par.nstep)
        self.field.f_xyz = np.zeros([discrete_t_par.nstep, 3])
        # Pi pulse: fmax*(cos^2(pi(t-t0)/(2s)) * cos(w(t-t0)))
        for i in range(discrete_t_par.nstep):
            if abs(discrete_t_par.dt * i - self.par.t0) < self.par.sigma:
                self.field.f_xyz[i] = self.par.fi \
                                      * np.square(
                    np.cos(np.pi * (discrete_t_par.dt * i - self.par.t0) / (2 * self.par.sigma))) \
                                      * np.cos(self.par.omega * (discrete_t_par.dt * i - self.par.t0))


    def sin_pulse(self, discrete_t_par):
        self.field.time_axis = np.linspace(0, discrete_t_par.dt * (discrete_t_par.nstep - 1), discrete_t_par.nstep)
        self.field.f_xyz = np.zeros([discrete_t_par.nstep, 3])
        # fmax*exp[-(t-t0)^2/s^2]*sin(wt)
        for i in range(discrete_t_par.nstep):
            self.field.f_xyz[i] = self.par.fi \
                                  * np.exp(-(discrete_t_par.dt * i - self.par.t0) ** 2 / (2 * (self.par.sigma) ** 2)) \
                                  * np.sin(self.par.omega * (discrete_t_par.dt * i))


    def gau_pulse(self, discrete_t_par):
        self.field.time_axis = np.linspace(0, discrete_t_par.dt * (discrete_t_par.nstep - 1), discrete_t_par.nstep)
        self.field.f_xyz = np.zeros([discrete_t_par.nstep, 3])
        for i in range(discrete_t_par.nstep):
            self.field.f_xyz[i] = self.par.fi \
                                  * np.exp(-(discrete_t_par.dt * i - self.par.t0) ** 2 / (self.par.sigma) ** 2)

    def sum_pulse(self, discrete_t_par):
        self.field.time_axis = np.linspace(0, discrete_t_par.dt * (discrete_t_par.nstep - 1), discrete_t_par.nstep)
        self.field.f_xyz = np.zeros([discrete_t_par.nstep, 3])
        # f0+fi*sin(wi*t)
        if self.par.omega.ndim == 1: #only one value of omega, 3D field [wx,wy,wz]
            for i in range(discrete_t_par.nstep):
                self.field.f_xyz[i] = self.par.fi[0] - self.par.fi[0] * np.sin(self.par.omega * i * discrete_t_par.dt) \
                                      + self.par.fi * np.sin(self.par.omega * i * discrete_t_par.dt)
        elif self.par.omega.ndim == 2: #matrix, N values of omega, 3D [w0x,w0y,w0z],[w1x,w1y,w1z]... il primo valore di omega deve essere sempre 0
            for i in range(discrete_t_par.nstep):
                self.field.f_xyz[i] = self.par.fi[0]
                for j in range(self.par.omega.shape[0]):
                    self.field.f_xyz[i] = self.field.f_xyz[i] \
                                          + self.par.fi[j] * np.sin(self.par.omega[j] * i * discrete_t_par.dt)



    def sin_cos_pulse(self, discrete_t_par):
        self.field.time_axis = np.linspace(0, discrete_t_par.dt * (discrete_t_par.nstep - 1), discrete_t_par.nstep)
        self.field.f_xyz = np.zeros([discrete_t_par.nstep, 3])
        # f0+fi*sin(wi*t)
        for i in range(discrete_t_par.nstep):
            self.field.f_xyz[i] = self.par.fi[0]
            for j in range(self.par.omega.shape[0]):
                    self.field.f_xyz[i] = self.field.f_xyz[i] + self.par.fi_cos[j] * np.cos(self.par.omega[j] * i * discrete_t_par.dt) - self.par.fi[j] * np.sin(self.par.omega[j] * i * discrete_t_par.dt)



    def genetic_pulse(self, discrete_t_par):
        self.chose_omega_fourier(discrete_t_par)
        self.par.fi = np.full([self.par.omega.shape[0], 3], self.par.fi[0])
        self.par.fi_cos = np.full([self.par.omega.shape[0], 3], self.par.fi[0])
        self.sum_pulse(discrete_t_par)


    def restart_genetic(self, discrete_t_par):
        self.sum_pulse(discrete_t_par)


    def chose_omega_fourier(self, discrete_t_par, n_fourier_omega = 0):
        if n_fourier_omega == 0:
            n_fourier_omega = int(2 + self.par.omega_sys[-1]*discrete_t_par.dt*discrete_t_par.nstep/(2*np.pi))
        self.par.omega = np.zeros([n_fourier_omega, 3])
        for n in range(n_fourier_omega):
            omega_n = (np.pi*2*n)/(discrete_t_par.dt*discrete_t_par.nstep)
            self.par.omega[n] = [omega_n, omega_n, omega_n]


    def chose_omega_energy(self):
        freq= []
        for i in range(self.par.omega_sys.shape[0]):
            for j in range(i):
                freq.append(self.par.omega_sys[i]-self.par.omega_sys[j])
        self.par.omega = np.zeros([len(freq), 3])
        self.par.omega[:,0] = freq
        self.par.omega[:, 1] = freq
        self.par.omega[:,2] = freq
