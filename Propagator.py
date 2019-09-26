import numpy as np
import abc

import auxiliary_functions as af
from Molecule import Molecule
from PCM import PCM, DinamicPCM, FrozenSolventPCM
import math_functions as mf


class Propagator():
    def __init__(self):
        self.propagator_terms = PropagatorTerms()
        self.propagator = []

    def set_propagator(self, dt, molecule, env):
        return

    def propagate_one_step(self, *args):
        return

    def propagate_n_step(self, *args):
        return

    def init_propagaror_terms(self, dt, molecule, pcm):
        self.propagator_terms.set_attributes(dt, molecule, pcm)
        self.propagator_terms.init_terms_dictionary()

    def add_term_to_propagator(self, term_name):
            self.propagator.append(self.propagator_terms.dict_terms[term_name])

    def clean_propagator(self):
        self.propagator = []





class PropagatorTerms():
    def __init__(self):
        self.dt = None
        self.mol = Molecule()
        self.pcm = PCM()

        self.dict_terms = {}

    def set_attributes(self, dt, molecule, pcm):
        self.dt = dt
        self.mol = molecule
        self.pcm = pcm

    def init_terms_dictionary(self):
        self.dict_terms["eulero1_coeff"] = self.eulero1_coeff_term
        self.dict_terms["eulero2_coeff"] = self.eulero1_coeff_term
        self.dict_terms["eulero_energy"] = self.eulero_energy_term
        self.dict_terms["eulero_field"] = self.eulero_field_term
        self.dict_terms["norm"] = self.norm
        self.dict_terms["oc_pcm_bwd"] = self.bwd_PCM_term
        self.dict_terms["eulero_pcm"] = self.eulero_Env_term

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

    def eulero_energy_term(self, order, *args):
        self.mol.wf.ci += -order * 1j * self.dt * (self.mol.en_ci * self.mol.wf.ci_prev[0])

    def eulero_field_term(self, order, field, *args):
        self.mol.wf.ci += -order * 1j * self.dt * (-np.dot(np.dot(self.mol.wf.ci_prev[0], self.mol.muT), field))

    def eulero_Env_term(self, order, field, *args):
        self.pcm.propagate(self.mol, field)
        self.mol.wf.ci += -order * 1j * self.dt \
                          * (np.dot(self.mol.wf.ci_prev[0],
                                    af.single_summation_tessere(self.pcm.get_q_t() - self.pcm.q00n, self.mol.Vijn)))

    def bwd_PCM_term(self, order, field, wf_fwd, *args):
        self.pcm.propagate_bwd_oc(self.mol, field, wf_fwd)
        self.mol.wf.ci += -order * 1j * self.dt \
                          * (np.dot(self.mol.wf.ci_prev[0],
                                    af.single_summation_tessere(self.pcm.get_q_t() - self.pcm.q00n,
                                                                self.mol.Vijn))
                             + np.dot(wf_fwd,
                                      af.single_summation_tessere(
                                            af.double_summation(self.mol.wf.ci_prev[0],
                                                                np.conj(wf_fwd),
                                                                self.mol.Vijn),
                                            self.pcm.qijn))
                             - np.dot(wf_fwd,
                                      af.single_summation_tessere(
                                            af.double_summation(wf_fwd,
                                                                np.conj(self.mol.wf.ci_prev[0]),
                                                                self.mol.Vijn),
                                            self.pcm.qijn)))




#    def bwd_PCM_term_fortran(self, order, field, wf_fwd, *args):
#        self.mol.wf.ci += -order * 1j * self.dt \
#                          * (mf.bwd_pcm(np.asfortranarray(self.mol.wf.ci_prev[0], dtype=np.complex128),
#                                        np.asfortranarray(wf_fwd, dtype=np.complex128),
#                                        np.asfortranarray(self.env.qijn_flip, dtype=np.complex128),
#                                        np.asfortranarray(self.env.Vij_flip, dtype=np.complex128),
#                                        np.asfortranarray(self.h.toSubtract, dtype=np.complex128)))







    #   def eulero_Env_TDPLas_term(self, order, field, *args):
    #       self.wavef.ci += -order*1j*self.dt \
    #                        * (sub_fortran.fwd_pcm(np.asfortranarray(self.wavef.ci_prev[0], dtype=np.complex128),
    #                                              np.asfortranarray(self.env.Vij_flip, dtype=np.complex128),
    #                                              np.asfortranarray(np.dot(field, self.env.cavity[:,0:3]), dtype=np.complex128),
    #                                              np.asfortranarray(self.env.to_subtract_fromH, dtype=np.complex128)))



#    def eulero_PCM_term_fortran(self, order, *args):
#        self.mol.wf.ci += -order * 1j * self.dt \
#                          * (mf.fwd_pcm(np.asfortranarray(self.mol.wf.ci_prev[0], dtype=np.complex128),
#                                        np.asfortranarray(self.env.qijn_flip, dtype=np.complex128),
#                                        np.asfortranarray(self.env.Vij_flip, dtype=np.complex128),
#                                        np.asfortranarray(self.env.to_subtract_fromH,dtype=np.complex128)))


#    def eulero_NANOP_term(self, order, *args):

#        daDareATDPlas=af.double_summation(self.wavef.ci_prev[0],
#                            np.conj(self.wavef.ci_prev[0]),
#                            self.env.get_Vijn())

#        qijnt_tdplas=cosa_presa_da_tdplas

#        self.wavef.ci += -order*1j*self.dt \
#                         * (np.dot(self.wavef.ci_prev[0],
#                                   af.single_summation_tessere(qijnt_tdplas, self.env.get_Vij())
#                                   -self.h.toSubtract))






