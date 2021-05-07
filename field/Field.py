import numpy as np

from SystemObj import DiscreteTimePar, Func_tMatrix

#---------------------------------------------
#field matrix at different times to be read or created
#---------------------------------------------
from parameters.FieldParameters import FieldParameters



#given the field parameters (included dt and nstep) and the field type fill self.field with the values of the field at each time step
#can read the field from file value by value, "read" choice, or read amplitudes and comegas and create a new field "sum" type

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
        self.par.nstep = field_input.nstep
        #additional steps where the field is equal to 0
        self.par.additional_steps = field_input.additional_steps

        dt = field_input.dt

        field = Func_tMatrix()
        field.f_xyz = field_input.field
        field.time_axis = field_input.field_time_axis
        self.chose_field(self.par.field_type, dt, field)



    def chose_field(self, key, dt = None, field = None):
         field_choice = {
            'const': lambda : self.const_pulse(dt),
            'pip': lambda : self.pip_pulse(dt),
            'sin': lambda: self.sin_pulse(dt),
            'gau': lambda: self.gau_pulse(dt),
            'sum': lambda: self.sum_pulse(dt),
            'genetic' : lambda: self.genetic_pulse(dt),
            'read': lambda : self.read_pulse(dt),
            'read_genetic': lambda : self.read_genetic_pulse(dt),
             #only internal values
            'restart_rabitz' : lambda: self.restart_rabitz(field),
            'restart_genetic' : lambda: self.restart_genetic(dt)
             }
         return field_choice.get(key, lambda: "Inexistent field type")()



    def add_empty_steps(self, dt):
        additional_t=np.linspace(dt * (self.par.nstep),
                                 dt * (self.par.nstep + self.par.additional_steps),
                                 self.par.additional_steps)
        additional_field = np.zeros((additional_t.shape[0],3))
        self.field.time_axis = np.concatenate((self.field.time_axis,additional_t))
        self.field.f_xyz = np.concatenate((self.field.f_xyz,additional_field))


    def restart_rabitz(self, field):
        self.field = field

    def read_pulse(self, dt):
        self.field.time_axis = np.loadtxt(self.par.restart_name, usecols=(2))
        print(self.field.time_axis.shape)
        self.field.f_xyz = np.loadtxt(self.par.restart_name, usecols = (3,4,5))


    def read_genetic_pulse(self, dt):
        load = np.loadtxt(self.par.restart_name) # w0, w1, ...wn \n a0 a1 ...an
        n = int((load.shape[0])/2)
        self.par.fi = load[n:]
        self.par.omega = load[:n]
        self.sum_pulse(dt)


    def const_pulse(self, dt):
        self.field.time_axis = np.linspace(0, dt * (self.par.nstep-1), self.par.nstep)
        self.field.f_xyz = np.zeros([self.par.nstep, 3])
        self.field.f_xyz[:, 0] = self.par.fi[0]
        self.field.f_xyz[:, 1] = self.par.fi[1]
        self.field.f_xyz[:, 2] = self.par.fi[2]
        if(self.par.additional_steps!=0):
            self.add_empty_steps(dt)


    def pip_pulse(self, dt):
        self.field.time_axis = np.linspace(0, dt * (self.par.nstep - 1), self.par.nstep)
        self.field.f_xyz = np.zeros([self.par.nstep, 3])
        # Pi pulse: fmax*(cos^2(pi(t-t0)/(2s)) * cos(w(t-t0)))
        for i in range(self.par.nstep):
            if abs(dt * i - self.par.t0) < self.par.sigma:
                self.field.f_xyz[i] = self.par.fi \
                                      * np.square(
                    np.cos(np.pi * (dt * i - self.par.t0) / (2 * self.par.sigma))) \
                                      * np.cos(self.par.omega * (dt * i - self.par.t0))
        if(self.par.additional_steps!=0):
            self.add_empty_steps(dt)

    def sin_pulse(self, dt):
        self.field.time_axis = np.linspace(0, dt * (self.par.nstep - 1), self.par.nstep)
        self.field.f_xyz = np.zeros([self.par.nstep, 3])
        # fmax*exp[-(t-t0)^2/s^2]*sin(wt)
        for i in range(self.par.nstep):
            self.field.f_xyz[i] = self.par.fi \
                                  * np.exp(-(dt * i - self.par.t0) ** 2 / (2 * (self.par.sigma) ** 2)) \
                                  * np.sin(self.par.omega * (dt * i))
        if(self.par.additional_steps!=0):
            self.add_empty_steps(dt)

    def gau_pulse(self, dt):
        self.field.time_axis = np.linspace(0, dt * (self.par.nstep - 1), self.par.nstep)
        self.field.f_xyz = np.zeros([self.par.nstep, 3])
        for i in range(self.par.nstep):
            self.field.f_xyz[i] = self.par.fi \
                                  * np.exp(-(dt * i - self.par.t0) ** 2 / (self.par.sigma) ** 2)
        if(self.par.additional_steps!=0):
            self.add_empty_steps(dt)


    def sum_pulse(self, dt):
        self.field.time_axis = np.linspace(0, dt * (self.par.nstep - 1), self.par.nstep)
        self.field.f_xyz = np.zeros([self.par.nstep, 3])
        # f0+fi*sin(wi*t)
        if self.par.omega.ndim == 1: #only one value of omega, 3D field [wx,wy,wz]
            for i in range(self.par.nstep):
                self.field.f_xyz[i] = self.par.fi[0] - self.par.fi[0] * np.sin(self.par.omega * i * dt) \
                                      + self.par.fi * np.sin(self.par.omega * i * dt)
        elif self.par.omega.ndim == 2: #matrix, N values of omega, 3D [w0x,w0y,w0z],[w1x,w1y,w1z]... il primo valore di omega deve essere sempre 0
            for i in range(self.par.nstep):
                self.field.f_xyz[i] = self.par.fi[0]
                for j in range(self.par.omega.shape[0]):
                    self.field.f_xyz[i] = self.field.f_xyz[i] \
                                          + self.par.fi[j] * np.sin(self.par.omega[j] * i * dt)
        if(self.par.additional_steps!=0):
            self.add_empty_steps(dt)



    def genetic_pulse(self, dt):
        self.chose_omega_fourier(dt)
        self.par.fi = np.full([self.par.omega.shape[0], 3], self.par.fi[0])
        self.sum_pulse(dt)


    def restart_genetic(self, dt):
        self.sum_pulse(dt)


    def chose_omega_fourier(self, dt, n_fourier_omega = 0):
        if n_fourier_omega == 0:
            n_fourier_omega = int(2 + self.par.omega_sys[-1]*dt*self.par.nstep/(2*np.pi))
        self.par.omega = np.zeros([n_fourier_omega, 3])
        for n in range(n_fourier_omega):
            omega_n = (np.pi*2*n)/(dt*self.par.nstep)
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
