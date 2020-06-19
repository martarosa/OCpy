# -*- coding: utf-8 -*-
"""
Created on Thu Mar 19 22:21:22 2020

@author: castd
"""
import numpy as np
import matplotlib.pyplot as plt
from field.Field import Field
import read_and_set.read.auxiliary_functions as af
from SystemObj import DiscreteTimePar
from astropy.convolution import Gaussian1DKernel, convolve
from mpl_toolkits import mplot3d

### import results ###

results = np.load("/Users/castd/Desktop/OCpy_da_spostare/test_scipy_quantum_10levresult.npy", allow_pickle = True)
fields_to_work_out = []
time_fields_to_work_out = []
field_fluency = []
iterations = 200

g = Gaussian1DKernel(stddev=1.5)



for k in range(iterations):
    fields_to_work_out.append(results[0].allvecs[5*(k)])
    
cost_function = np.asarray(results[1][1::5])

#####

### create field and generate fourier_transform and save ###


n_points_field = 250
dt = 1
discrete_t_par = DiscreteTimePar()
discrete_t_par.dt = dt
discrete_t_par.nstep = n_points_field
field = Field()
field.par.omega_sys = np.append([0,], np.loadtxt("/Users/castd/Desktop/OCpy_da_spostare/ci_energy.inp", usecols=((3,)))/27.2114)
field.chose_omega_energy()
field.chose_omega_fourier(discrete_t_par=discrete_t_par)
print(field.par.omega_sys)

def alpha_field_J_integral(field):
    ax_square= field.field.f_xyz.ndim - 1
    ax_integral= field.field.f_xyz.ndim - 2
    f_square = np.sum(field.field.f_xyz * field.field.f_xyz, axis=ax_square) 
    f_integral = np.sum(f_square, axis=ax_integral)
    out_integral = f_integral*discrete_t_par.dt
    return out_integral

#nharmonic = len(omega_fourier)
N_points_zero = 500
N_tot = N_points_zero + n_points_field

iteration_freq_matrix = np.zeros((iterations, len(np.linspace(0.0, np.pi*2/(N_tot*dt)*N_tot, N_tot)), 2))
#iteration_freq_matrix = np.zeros((iterations, len(freq_x_axis)))
#z_iteration_freq_matrix = np.zeros((iterations, freq_x_axis))

for k in range(iterations):
    amplitudes = fields_to_work_out[k]
    field.par.fi = np.asarray(amplitudes).reshape((-1, 3))
    print(field.par.fi)
    field.chose_field('sum', discrete_t_par = discrete_t_par)
    field_fluency.append(alpha_field_J_integral(field))
    for j in range(2):
        iteration_freq_matrix[k,:,j] = convolve(af.smooth_fft_field(dt, n_points_field, N_points_zero, field.field.f_xyz[:,j])[1][:], g)
        
pop_tgt = 1 + np.asarray(field_fluency) - cost_function

new_cost_function = pop_tgt - np.asarray(field_fluency)

#
#plt.rcParams['figure.figsize'] = 17.5, 5.0
#plt.rc('font', **{'family': 'sans-serif', 'sans-serif': ['Computer Modern'], 'size' : 16})
#plt.rc('lines', linewidth=1.5)
#
#fig = plt.figure()
#
#ax1 = fig.add_subplot(131)
#
#ax1.plot(np.arange(0,1000,5), new_cost_function, c = 'C3')
#
#ax1.set_xlim(0,1000)
#ax1.set_ylim(0,1)
#ax1.set_yticks(np.arange(0,1.2,0.2))
#
#ax1.set_xlabel("Iterations")
#ax1.set_ylabel("$J[a]$")
#
#    
#ax2 = fig.add_subplot(132)
#
#ax2.plot(np.arange(0,1000,5), pop_tgt, c = 'C3')
#
#ax2.set_xlim(0,1000)
#ax2.set_ylim(0,1)
#ax2.set_yticks(np.arange(0,1.2,0.2))
#
#ax2.set_xlabel("Iterations")
#ax2.set_ylabel("$P_{target}(T)$")
#
#ax3 = fig.add_subplot(133)
#
#ax3.plot(np.arange(0,1000,5), field_fluency, c = 'C3')
#
#ax3.set_xlim(0,1000)
#ax3.set_ylim(0,0.2)
#ax3.set_yticks(np.arange(0,0.25,0.05))
#
#ax3.set_xlabel("Iterations")
#ax3.set_ylabel("$\int_0^T\alpha(t)|E(t)|^2dt$")
#
#plt.tight_layout()
#plt.savefig('/Users/castd/Desktop/cost_function', dpi=600)
#plt.show()

####

x_iteration_freq_matrix = iteration_freq_matrix[:,:,0]
y_iteration_freq_matrix = iteration_freq_matrix[:,:,1]


#surface plot of populations as function of sites and time
fig = plt.figure()
ax = plt.axes(projection='3d')
X, Y = np.meshgrid(np.linspace(0.0, np.pi*2/(N_tot*dt)*N_tot, N_tot)[:80], np.arange(iterations))
ax.plot_surface(X, Y, y_iteration_freq_matrix[:,:80], rstride=1, cstride=1, cmap='viridis', edgecolor='none')
ax.view_init(azim=85, elev=30)
ax.set_ylabel(r'Iterations', labelpad = 15)
ax.set_xlabel(r'$\omega [a.u.] $', labelpad = 20)
ax.set_zlabel(r'$\mathcal{F}(\omega)$')
ax.set_zticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
ax.set_ylim3d(0,200)
ax.set_yticks([int(element)-1 for element in np.linspace(0,200,1000)[5::200]])
ax.set_xlim3d(0,np.linspace(0.0, np.pi*2/(N_tot*dt)*N_tot, N_tot)[:80][-1])
#ax.invert_xaxis()
#plt.tight_layout()
plt.savefig('/Users/castd/Desktop/3dplot_trial', dpi=600)
plt.show()




