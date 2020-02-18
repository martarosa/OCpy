import numpy as np

from propagator import math_functions as mf

from pcm.ABCPCM import ABCPCM, PCMParameters
from read_and_set.read import auxiliary_functions as af

# if the solvent is frozen qijn and qijn_lf are fixed, calculated once at the initialization
#q_t dependence on time comes only from <psi(t)|qijn|psi(t)>
#q_t_lf dependence comes from qijn_lf * field_t
class FrozenSolventPCM(ABCPCM):
    def __init__(self):
        super().__init__()
        self.par = PCMParameters()

        self.q_t = None
        self.qijn = None
        self.qijn_lf = None
        self.q00n = None

        self.qijn_fortran_flip = None


    def init_pcm(self, PCM_input, mol, field_dt_vector):
        self.par.env = PCM_input.env
        self.par.cavity = PCM_input.cavity
        self.init_static_matrices(PCM_input.Qnn_reactionfield, PCM_input.Qnn_localfield, mol)
        self.par.muLF = -af.matrix_prod_tesserae_ijn_nn(self.qijn_lf, mol.par.Vijn)


        self.propagate(0, mol, field_dt_vector)
        self.qijn_fortran_flip = af.flip_3D_py2f(self.qijn)


    def propagate(self, i, mol, field_dt_vector):
        q_t_reactionf = af.double_summation(mol.wf.ci_prev[0], np.conj(mol.wf.ci_prev[0]), self.qijn)
        q_t_lf = np.dot(self.qijn_lf, field_dt_vector)
        self.q_t = np.array([q_t_reactionf, q_t_lf])


    def propagate_fortran(self, i, mol, field_dt_vector):
        q_t_reactionf = mf.propagate_q_frozen(mol.wf.ci_prev[0], self.qijn_fortran_flip)
        q_t_lf = np.dot(self.qijn_lf, field_dt_vector)
        self.q_t = np.array([q_t_reactionf, q_t_lf])


    def propagate_bwd_oc(self, i, chi_ci, field_dt_vector):
        q_t_reactionf = af.double_summation(chi_ci, np.conj(chi_ci), self.qijn)
        q_t_lf = np.dot(self.qijn_lf, field_dt_vector)
        self.q_t = np.array([q_t_reactionf, q_t_lf])

    def propagate_bwd_oc_fortran(self, i, chi_ci, field_dt_vector):
        q_t_reactionf = mf.propagate_q_frozen(chi_ci, self.qijn_fortran_flip)
        q_t_lf = np.dot(self.qijn_lf, field_dt_vector)
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
        self.calc_qijn(Q_gamess, mol.par.Vijn)
        self.calc_qijn_lf(Q_gamess_lf)

    def calc_qijn(self, Q_gamess, Vijn):
        Q_tdplas = np.dot(Q_gamess, np.diag(self.par.cavity[:,3]))
        self.qijn = af.matrix_prod_tesserae_ijn_nn(Q_tdplas, Vijn)
        self.q00n = self.qijn[0,0]

    def calc_qijn_lf(self, Q_gamess_lf):
        Q_local_field_tdplas = np.dot(Q_gamess_lf, np.diag(self.par.cavity[:,3]))
        qijn_lf = - np.dot(Q_local_field_tdplas, self.par.cavity[:,0:3])
        self.qijn_lf = qijn_lf.astype(float)