import numpy as np

from propagator import math_functions as mf

from medium.ABCMedium import ABCMedium
from parameters.MediumParameters import MediumParameters
from read_and_set.read import auxiliary_functions as af

#This approximation can be used only for the SOLVENT, not for the NP
# When the solvent is approximated as frozen, qijn and qijn_lf are fixed in time, calculated once when the
# medium is initialized (init_static_matrices).
#the propagation is the product between the fixed charges matrices and the wf coefficients multiplied by the reaction field
#q_t dependence on time comes only from <psi(t)|qijn|psi(t)>
#and ny the external pulse for the local field
#q_t_lf dependence comes from qijn_lf * field_t
#matrices products can be done in fortran to speed up the calculations, but it is still in OCpy, TDPlas is not used

class FrozenSolventMedium(ABCMedium):
    def __init__(self):
        super().__init__()
        self.par = MediumParameters()
        self.muLF = None
        # static matrices in frozen solvent, initialized once
        self.qijn = None
        self.qijn_lf = None
        self.q00n = None
        self.qijn_fortran_flip = None


    def init_medium(self, medium_input, mol, field_object):
        self.par.medium = medium_input.medium
        self.par.cavity = medium_input.cavity
        self.init_static_matrices(medium_input.Qnn_reactionfield, medium_input.Qnn_localfield, mol)
        self.muLF = -af.matrix_prod_tesserae_ijn_nn(self.qijn_lf, mol.par.Vijn)
        self.propagate_charges(mol, field_object.f_xyz[0])
        self.qijn_fortran_flip = af.flip_3D_py2f(self.qijn)


    def reset_medium(self, *args):
        pass


    def propagate_charges(self, mol, field_dt_vector):
        q_t_reactionf = af.double_summation(mol.wf.ci_prev[0], np.conj(mol.wf.ci_prev[0]), self.qijn)
        q_t_lf = np.dot(self.qijn_lf, field_dt_vector)
        self.par.q_t = q_t_reactionf + q_t_lf


    def propagate_charges_fortran(self, mol, field_dt_vector):
        q_t_reactionf = mf.propagate_q_frozen(mol.wf.ci_prev[0], self.qijn_fortran_flip)
        q_t_lf = np.dot(self.qijn_lf, field_dt_vector)
        self.par.q_t = q_t_reactionf + q_t_lf


    def propagate_bwd_oc(self, chi_ci, field_dt_vector):
        q_t_reactionf = af.double_summation(chi_ci, np.conj(chi_ci), self.qijn)
        q_t_lf = np.dot(self.qijn_lf, field_dt_vector)
        self.par.q_t = q_t_reactionf + q_t_lf

    def propagate_bwd_oc_fortran(self, chi_ci, field_dt_vector):
        q_t_reactionf = mf.propagate_q_frozen(chi_ci, self.qijn_fortran_flip)
        q_t_lf = np.dot(self.qijn_lf, field_dt_vector)
        self.par.q_t = q_t_reactionf + q_t_lf


    def get_q_t(self):
        return self.par.q_t - self.q00n


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

