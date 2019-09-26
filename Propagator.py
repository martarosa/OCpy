import numpy as np
import abc

import auxiliary_functions as af
from Molecule import Molecule
from PCM import PCM, DinamicPCM, FrozenSolventPCM
import math_functions as mf

class Propagator():
    def __init__(self):
        self.dt = None
        self.mol = Molecule()
        self.pcm = PCM()

        self.propagator = []

    def set_propagator(self, dt, molecule, env):
        return

    def propagate_one_step(self, *args):
        return

    def propagate_n_step(self, *args):
        return


    def set_attributes(self, dt, molecule, pcm):
        self.dt = dt
        self.mol = molecule
        self.pcm = pcm

    def test(self):
        print("manamanamana")


    def clean_propagator(self):
        self.propagator = []

    def add_term_to_propagator(self, term_name):
        if term_name == "eulero1_coeff":
            self.propagator.append(self.eulero1_coeff_term)
        if term_name == "eulero2_coeff":
            self.propagator.append(self.eulero2_coeff_term)
        if term_name == "eulero_energy":
            self.propagator.append(self.eulero_energy_term)
        elif term_name == "eulero_field":
            self.propagator.append(self.eulero_field_term)
        elif term_name == "norm":
            self.propagator.append(self.mol.wf.norm_ci)



    # <editor-fold desc="H terms">
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

    #   def eulero_Env_TDPLas_term(self, order, field, *args):
    #       self.wavef.ci += -order*1j*self.dt \
    #                        * (sub_fortran.fwd_pcm(np.asfortranarray(self.wavef.ci_prev[0], dtype=np.complex128),
    #                                              np.asfortranarray(self.env.Vij_flip, dtype=np.complex128),
    #                                              np.asfortranarray(np.dot(field, self.env.cavity[:,0:3]), dtype=np.complex128),
    #                                              np.asfortranarray(self.env.to_subtract_fromH, dtype=np.complex128)))

    def eulero_Env_term(self, order, field, *args):
        self.pcm.propagate(self.mol, field)
        self.mol.wf.ci += -order * 1j * self.dt \
                          * (np.dot(self.mol.wf.ci_prev[0],
                                    af.single_summation_tessere(self.pcm.get_q_t() - self.pcm.q00n, self.mol.Vijn)))

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


