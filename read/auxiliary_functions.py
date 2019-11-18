import numpy as np
from scipy.fftpack import fft
import os


def normalize_vector(vector):
    normalized = vector/np.sqrt(np.dot(np.conj(vector), vector))
    return normalized

    #|t><t|v>
def apply_projection(vector, vector_to_project):
    projected = np.dot(np.outer(vector_to_project, np.conj(vector_to_project)), vector)
    return projected

    #<v|t><t|v> = <v|O|v>
def projector_mean_value(vector, vector_to_project):
    mean_value = np.dot(np.conj(vector), apply_projection(vector, vector_to_project))
    return mean_value


def normalized_projector_mean_value(vector, vector_to_project):
    mean_value = np.dot(np.conj(vector), apply_projection(vector, vector_to_project))/np.dot(vector,vector)
    return mean_value


def field_J_integral(field_matrix, dt):
    ax_square=field_matrix.ndim-1
    ax_integral=field_matrix.ndim-2
    f_square = np.sum(field_matrix*field_matrix,axis=ax_square)
    f_integral = np.sum(f_square, axis=ax_integral)
    out_field = f_integral*dt
    return out_field


def alpha_field_J_integral(field_matrix, alpha_t, dt):
    ax_square=field_matrix.ndim-1
    ax_integral=field_matrix.ndim-2
    f_square = np.sum(field_matrix*field_matrix,axis=ax_square)*alpha_t
    f_integral = np.sum(f_square, axis=ax_integral)
    out_integral = f_integral*dt
    return out_integral


def smooth_fft_field(time_step, N_points, N_points_zero, field):
    zero = np.zeros([N_points_zero])
    longer_field = np.hstack((field, zero))
    N_tot = N_points_zero + N_points
    xf = np.linspace(0.0, np.pi*2/(N_tot*time_step)*N_tot, N_tot)
    yf = fft(longer_field)
    fourier_plot = np.array([xf, np.abs(yf)])
    return fourier_plot



def single_summation_tessere(ket, M):
    single_sum = np.sum(M*ket[np.newaxis, np.newaxis, :], axis=2)
    return single_sum


def double_summation(ket, bra, M):
    double_sum = np.dot(bra, np.sum(M*ket[np.newaxis, :, np.newaxis], axis=1))
    return double_sum


def double_summation_state_tessere(Mijn,Nijn):
    tmp=[]
    for i in range(Mijn.shape[2]):
        tmp.append(np.dot(Mijn[:,:,i],Nijn[:,:,i]))
    tmp=np.asarray(tmp)
    out=np.sum(tmp,axis=0)
    return out


def matrix_prod_tesserae_ijn_ijn(Mijn,Nijn):
    states=Mijn.shape[0]
    tessere=Mijn.shape[2]
    out=np.zeros((states,states,states))
    for i in range(states):
        for j in range(states):
            tmp=[]
            for n in range(tessere):
                tmp.append(Mijn[i,:,n]*Nijn[:,j,n])
            tmp=np.asarray(tmp)
            out[i,j] = tmp.sum(axis=0)
    return out

def matrix_prod_tesserae_ijn_nn(Mnm, Mijn):
    prod = np.dot(Mijn, Mnm)
    return prod.astype(float)



def func(Mnr,Mijn):
    prod = np.dot(Mijn, Mnr)
    return prod



def population_from_wf_vector(wf):
    out = (np.conj(wf)*wf)/np.dot(np.conj(wf), wf)
    return out


def population_from_wf_matrix(wf):
    out = np.apply_along_axis(population_from_wf_vector, 1, wf)
    return out




def remember_path(namefile):
    folder = os.path.split(str(namefile))[0]
    folder = folder.split("='")[1]
    return folder


def convert_pot_to_waveT_format(matrix, name_out_file, n_en):
    V_ijn_el = matrix[1:, 1:]
    il = np.tril_indices(V_ijn_el.shape[0])
    V_trilled = V_ijn_el[il]
    tril_index = 0
    for i in range(V_ijn_el.shape[0]):
        for j in range(V_ijn_el.shape[1]):
            if (j <= i):
                f = open(name_out_file, 'a')
                f.write(str(i + 1) + ' ' + str(j + 1) + '\n')
                f.close()
                f = open(name_out_file, 'ab')
                np.savetxt(f, V_trilled[tril_index], fmt='%10.5e')
                tril_index += 1
                f.close()





def calc_deltaJ_analitical_terms(psi_k, psi_km1, chi_k, eps_k, eps_km1, eps_tilde_k, Q, V, alpha, dt, nstep, ntessere):
    t1 = alpha_field_J_integral(np.real(eps_k-eps_tilde_k), alpha, dt)
    t2 = alpha_field_J_integral(np.real(eps_tilde_k-eps_km1), alpha, dt)
    t3_t = np.zeros((nstep, ntessere), dtype=complex)
    t4_t = np.zeros((nstep, ntessere), dtype=complex)
    for i in range(nstep):
        t3_t[i] = (2*np.imag(double_summation(psi_k[i], np.conj(chi_k[i]), V))*(double_summation(psi_k[i], np.conj(psi_k[i]), Q)-
                                                              double_summation(psi_km1[i],np.conj(psi_km1[i]),Q)))
        t4_t[i] = (2*np.imag(double_summation(psi_km1[i], np.conj(psi_k[i]), Q))
               *2*np.real(double_summation(psi_km1[i], np.conj(chi_k[i]), V)))
    t3 = np.sum(t3_t)*dt
    t4 = np.sum(t4_t)*dt
    deltaJ = t1+t2+t3+t4
    return deltaJ

def flip_3D_py2f(Mpy):
    tmp = np.swapaxes(Mpy, 1, 2)
    Mf = np.swapaxes(tmp, 0, 1)
    return Mf


def flip_3D_f2py(Mf):
    tmp = np.swapaxes(Mf, 0, 1)
    Mpy = np.swapaxes(tmp, 1, 2)
    return Mpy
