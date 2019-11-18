import numpy as np
from read import auxiliary_functions as af
from propagator import math_functions as mf

class PCM():
    def __init__(self):
        self.env = None
        self.cavity = None
        self.q_t = None

        self.muLF = 0
        self.qijn = None
        self.qijn_lf = None
        self.q00n = None

    def init_pcm(self, PCM_parameters, mol, field_t):
        pass

    def propagate(self, mol, field_t):
        pass

    def get_q_t(self):
        pass



# if the solvent is frozen qijn and qijn_lf are fixed, calculated once at the initialization
#q_t dependence on time comes only from <psi(t)|qijn|psi(t)>
#q_t_lf dependence comes from qijn_lf * field_t
class FrozenSolventPCM(PCM):
    def __init__(self):
        super().__init__()
        self.env = None
        self.cavity = None
        self.q_t = None

        self.muLF = 0
        self.qijn = None
        self.qijn_lf = None
        self.q00n = None

        self.qijn_fortran_flip = None
        self.qijn_lf_fortran_flipped = None


    def init_pcm(self, PCM_parameters, mol, field_t):
        self.env = PCM_parameters.env
        self.cavity = PCM_parameters.cavity
        self.init_static_matrices(PCM_parameters.Qnn_reactionfield, PCM_parameters.Qnn_localfield, mol)
        self.muLF = -af.matrix_prod_tesserae_ijn_nn(self.qijn_lf, mol.Vijn)
        self.propagate(mol, field_t)
        self.qijn_fortran_flip = af.flip_3D_py2f(self.qijn)

    def propagate(self, mol, field_t):
        q_t_reactionf = af.double_summation(mol.wf.ci_prev[0], np.conj(mol.wf.ci_prev[0]), self.qijn)
        q_t_lf = np.dot(self.qijn_lf, field_t)
        self.q_t = np.array([q_t_reactionf, q_t_lf])

    def propagate_fortran(self, ci, field_t):
        q_t_reactionf = mf.propagate_q_frozen(ci, self.qijn_fortran_flip)
        q_t_lf = np.dot(self.qijn_lf, field_t)
        self.q_t = np.array([q_t_reactionf, q_t_lf])


    def propagate_bwd_oc(self, mol, field_t, ci):
        q_t_reactionf = af.double_summation(ci, np.conj(ci), self.qijn)
        q_t_lf = np.dot(self.qijn_lf, field_t)
        self.q_t = np.array([q_t_reactionf, q_t_lf])

    def get_q_t(self):
        # q_tmp= self.q_t.q_t[0]+self.q_t.q_t[1]
        # diff= 0 - q_tmp
        # print(diff)
        # delta_x_tessera = diff/self.q00n.shape[-1]
        # q_tot = q_tmp + delta_x_tessera
        # return q_tot
        return self.q_t[0] + self.q_t[1]

    def init_static_matrices(self, Q_gamess, Q_gamess_lf, mol):
        self.calc_qijn(Q_gamess, mol.Vijn)
        self.calc_qijn_lf(Q_gamess_lf)

    def calc_qijn(self, Q_gamess, Vijn):
        Q_tdplas = np.dot(Q_gamess, np.diag(self.cavity[:,3]))
        self.qijn = af.matrix_prod_tesserae_ijn_nn(Q_tdplas, Vijn)
        self.q00n = self.qijn[0,0]

    def calc_qijn_lf(self, Q_gamess_lf):
        Q_local_field_tdplas = np.dot(Q_gamess_lf, np.diag(self.cavity[:,3]))
        qijn_lf = - np.dot(Q_local_field_tdplas, self.cavity[:,0:3])
        self.qijn_lf = qijn_lf.astype(float)




# if the solvent is dinamic qijn and qijn_lf depend on time = qijn_t, qijn_lf_t
#q_t dependence on time comes from <psi(t)|qijn_t|psi(t)>
#q_t_lf dependence comes from qijn_lf_t * field_t
#at each timestep qijn_t and qijn_lf_t are given by tdplas
class DinamicPCM(PCM):
    def __init__(self):
        super().__init__()
        self.env = None
        self.cavity = None
        self.q_t = None

    def init_pcm(self, PCM_parameters, mol, field_t):
        self.cavity = PCM_parameters.cavity
        input_V_tdplas = af.double_summation(mol.wf.ci_prev[0], np.conj(mol.wf.ci_prev[0]), mol.Vijn)
        input_V_lf_tdplas = - np.dot(field_t, self.cavity[:,0:3])
        interfaccia_tdplas.init_mdm()


    def propagate(self, mol, field_t):
        input_V_tdplas = af.double_summation(mol.wf.ci_prev[0], np.conj(mol.wf.ci_prev[0]), mol.Vijn)
        input_V_lf_tdplas = - np.dot(field_t, self.cavity[:,0:3])
        interfaccia_tdplas.propagate_mdm()


    def get_q_t(self):
        # q_tmp= self.q_t.q_t[0]+self.q_t.q_t[1]
        # diff= 0 - q_tmp
        # print(diff)
        # delta_x_tessera = diff/self.q00n.shape[-1]
        # q_tot = q_tmp + delta_x_tessera
        # return q_tot
        return self.q_t[0] + self.q_t[1]