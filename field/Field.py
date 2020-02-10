import numpy as np





#---------------------------------------------
#field matrix at different times to be read or created
#---------------------------------------------

class FieldParameters():
    def __init__(self):
        self.field_type = None
        self.fi = None
        self.omega = None
        self.sigma = None
        self.t0 =  None
        self.omega_sys = None


class Field():
    def __init__(self):
        self.field = None
        self.field_time_axis = None
        self.par = FieldParameters()


    def init_field(self, field_input):
        self.par.field_type = field_input.field_type
        self.par.fi = field_input.fi
        self.par.omega = field_input.omega
        self.par.sigma = field_input.sigma
        self.par.t0 = field_input.t0
        self.par.omega_sys = field_input.omega_sys

        self.chose_field(self.par.field_type,
                         field_input.dt,
                         field_input.nstep,
                         field_input.field,
                         field_input.field_time_axis)

    def chose_field(self, key, dt, nstep, field, field_time_axis):
         field = {
            'const': lambda : self.const_pulse(dt, nstep),
            'pip': lambda : self.pip_pulse(dt, nstep),
            'sin': lambda: self.sin_pulse(dt, nstep),
            'gau': lambda: self.gau_pulse(dt, nstep),
            'sum': lambda: self.sum_pulse(dt, nstep),
            'genetic' : lambda: self.genetic_pulse(dt, nstep),
             #only internal values
            'restart_rabitz' : lambda: self.restart_rabitz(field, field_time_axis),
            'restart_genetic' : lambda: self.sum_pulse(dt, nstep)
             }
         return field.get(key, lambda: "Inexistent field type")()



    def restart_rabitz(self, field, field_time_axis):
        self.field = field
        self.field_time_axis = field_time_axis



    def const_pulse(self, dt, nstep):
        self.field_time_axis = np.linspace(0, dt * nstep, nstep)
        self.field = np.zeros([nstep, 3])
        self.field[:, 0] = self.par.fi[0]
        self.field[:, 1] = self.par.fi[1]
        self.field[:, 2] = self.par.fi[2]




    def pip_pulse(self, dt, nstep):
        self.field_time_axis = np.linspace(0, dt * nstep, nstep)
        self.field = np.zeros([nstep, 3])
        # Pi pulse: fmax*(cos^2(pi(t-t0)/(2s)) * cos(w(t-t0)))
        for i in range(nstep):
            if abs(dt * i - self.par.t0) < self.par.sigma:
                self.field[i] = self.par.fi \
                                * np.square(
                    np.cos(np.pi * (dt * i - self.par.t0) / (2 * self.par.sigma))) \
                                * np.cos(self.par.omega * (dt * i - self.par.t0))


    def sin_pulse(self, dt, nstep):
        self.field_time_axis = np.linspace(0, dt * nstep, nstep)
        self.field = np.zeros([nstep, 3])
        # fmax*exp[-(t-t0)^2/s^2]*sin(wt)
        for i in range(nstep):
            self.field[i] = self.par.fi \
                            * np.exp(-(dt * i - self.par.t0) ** 2 / (2 * (self.par.sigma) ** 2)) \
                            * np.sin(self.par.omega * (dt * i))


    def gau_pulse(self, dt, nstep):
        self.field_time_axis = np.linspace(0, dt * nstep, nstep)
        self.field = np.zeros([nstep, 3])
        for i in range(nstep):
            self.field[i] = self.par.fi \
                            * np.exp(-(dt * i - self.par.t0) ** 2 / (self.par.sigma) ** 2)

    def sum_pulse(self, dt, nstep):
        self.field_time_axis = np.linspace(0, dt * nstep, nstep)
        self.field = np.zeros([nstep, 3])
        # f0+fi*sin(wi*t)
        if self.par.omega.ndim == 1: #only one value of omega, 3D field [wx,wy,wz]
            for i in range(nstep):
                self.field[i] = self.par.fi[0] - self.par.fi[0] * np.sin(self.par.omega * i * dt) \
                                + self.par.fi * np.sin(self.par.omega * i * dt)
        elif self.par.omega.ndim == 2: #matrix, N values of omega, 3D [w0x,w0y,w0z],[w1x,w1y,w1z]... il primo valore di omega deve essere sempre 0
            for i in range(nstep):
                self.field[i] = self.par.fi[0]
                for j in range(self.par.omega.shape[0]):
                    self.field[i] = self.field[i] \
                                    + self.par.fi[j] * np.sin(self.par.omega[j] * i * dt)


    def genetic_pulse(self,dt, nstep):
        self.field_time_axis = np.linspace(0, dt * nstep, nstep)
        self.chose_omega_fourier(dt, nstep)
        self.par.fi = np.full([self.par.omega.shape[0]+1, 3], self.par.fi[0])
        self.sum_pulse(dt, nstep)


    def chose_omega_fourier(self, dt, nstep):
        n_fourier_omega_min = 0
        n_fourier_omega_max = int(2 + self.par.omega_sys[-1]*dt*nstep/(2*np.pi))
        n_fourier_omega = n_fourier_omega_max - n_fourier_omega_min
        self.par.omega = np.zeros([n_fourier_omega, 3])
        for n in range(n_fourier_omega):
            omega_n = (np.pi*2*(n + n_fourier_omega_min))/(dt*nstep)
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

    def get_full_field(self):
        out_field = np.zeros([self.field.shape[0],4])
        out_field[:,0] = self.field_time_axis
        out_field[:,1:] = self.field
        return out_field