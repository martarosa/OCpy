import numpy as np

from read_and_set.read import auxiliary_functions as af
from molecule.Molecule import Molecule
from propagator import math_functions as mf

# all posible propagator term. Some use PCM methods specific of FrozenSolventPCM child class, If/when
# DinamicPCM(PCM) will have rabitz implemented, it will have a propagate_bwd_oc term, while the propagation in
# DinamicPCM(PCM) is always done with fortran, so propagate_fortran doesn't exist
# terms are added to the propagator that is used for the wave function propagation
class PropagatorTerms():
    def __init__(self):
        self.mol = Molecule()
        self.medium = None
        self.dict_terms = {}


    def set_attributes(self, molecule, medium):
        self.mol = molecule
        self.medium = medium

    def init_terms_dictionary(self):
        self.dict_terms["eulero1_coeff"] = self.eulero1_coeff_term
        self.dict_terms["eulero2_coeff"] = self.eulero2_coeff_term
        self.dict_terms["eulero_energy"] = self.eulero_energy_term
        self.dict_terms["eulero_field"] = self.eulero_field_term
        self.dict_terms["eulero_medium"] = self.eulero_medium_term_fortran
        self.dict_terms["norm"] = self.norm
        self.dict_terms["oc_medium_bwd"] = self.bwd_medium_term_fortran

    # <editor-fold desc="H terms">

    def norm(self, *args):
        self.mol.wf.norm_ci()

    def eulero1_coeff_term(self, *args):
        self.mol.wf.ci_prev[0] = np.copy(self.mol.wf.ci)
        self.mol.wf.ci = np.copy(self.mol.wf.ci_prev[0])

    def eulero2_coeff_term(self, *args):
        self.mol.wf.ci_prev[1] = np.copy(self.mol.wf.ci_prev[0])
        self.mol.wf.ci_prev[0] = np.copy(self.mol.wf.ci)
        self.mol.wf.ci = np.copy(self.mol.wf.ci_prev[1])

    def eulero_energy_term(self, order, dt, *args):
        self.mol.wf.ci += -order * 1j * dt * (self.mol.par.en_ci * self.mol.wf.ci_prev[0])

    def eulero_field_term(self, order, dt, field_dt_vector, *args):
        self.mol.wf.ci += -order * 1j * dt * (-np.dot(np.dot(self.mol.wf.ci_prev[0], self.mol.par.muT), field_dt_vector))


    def eulero_medium_term(self, order, dt, field_dt_vector, *args):
        self.medium.propagate(self.mol, field_dt_vector)
        self.mol.wf.ci += -order * 1j * dt \
                          * (np.dot(self.mol.wf.ci_prev[0],
                                    af.single_summation_tessere(self.medium.get_q_t() - self.medium.q00n, self.mol.par.Vijn)))


    def bwd_medium_term(self, order, dt, field_dt_vector, wf_fwd, *args):
        self.medium.propagate_bwd_oc(wf_fwd, field_dt_vector)
        self.mol.wf.ci += -order * 1j * dt \
                          * (np.dot(self.mol.wf.ci_prev[0],
                                    af.single_summation_tessere(self.medium.get_q_t() - self.medium.q00n,
                                                                self.mol.par.Vijn))
                             + np.dot(wf_fwd,
                                      af.single_summation_tessere(
                                            af.double_summation(self.mol.wf.ci_prev[0],
                                                                np.conj(wf_fwd),
                                                                self.mol.par.Vijn),
                                            self.medium.qijn))
                             - np.dot(wf_fwd,
                                      af.single_summation_tessere(
                                            af.double_summation(wf_fwd,
                                                                np.conj(self.mol.wf.ci_prev[0]),
                                                                self.mol.par.Vijn),
                                            self.medium.qijn)))


    def eulero_medium_term_fortran(self, order, dt, field_dt_vector, *args):
        self.medium.propagate_fortran(self.mol, field_dt_vector)
        q_t = self.medium.get_q_t()
        self.mol.wf.ci += -order * 1j * dt \
                          *(mf.eulero_pcm(np.asfortranarray(self.mol.wf.ci_prev[0], dtype=np.complex128),
                                         np.asfortranarray(q_t, dtype=np.complex128),
                                         np.asfortranarray(self.mol.par.Vijn_fortran_flip, dtype=np.complex128)))



    def bwd_medium_term_fortran(self, order, dt, field_dt_vector, wf_fwd, *args):
        self.medium.propagate_bwd_oc_fortran(wf_fwd, field_dt_vector)
        q_t = self.medium.get_q_t() - self.medium.q00n
        self.mol.wf.ci += -order * 1j * dt \
                          * (mf.bwd_pcm(np.asfortranarray(self.mol.wf.ci_prev[0], dtype=np.complex128),
                                        np.asfortranarray(wf_fwd, dtype=np.complex128),
                                        np.asfortranarray(q_t, dtype=np.complex128),
                                        np.asfortranarray(self.mol.par.Vijn_fortran_flip, dtype=np.complex128),
                                        np.asfortranarray(self.medium.qijn_fortran_flip, dtype=np.complex128)))



