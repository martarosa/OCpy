import numpy as np

from Propagator import Propagator
import auxiliary_functions as af

#import math_functions as mf

class PropagatorOCfwd(Propagator):
    def __init__(self):
        super().__init__()

    def set_propagator(self, dt, molecule, env):
        self.set_attributes(dt, molecule, env)
        self.clean_propagator()
        self.add_term_to_propagator("eulero1_coeff")
        self.add_term_to_propagator("eulero_energy")
        self.add_term_to_propagator("eulero_field")
        if self.pcm.env == "sol":
            self.add_term_to_propagator("eulero_pcm")

    def propagate_one_step(self, field):
        for func in self.propagator:
            func(1, field)

    def propagate_n_step(self, nstep, field):
        out = list()
        out.append(self.mol.wf.ci)
        for i in range(nstep):
            self.propagate_one_step(field[i])
            out.append(self.mol.wf.ci)
        out = np.array(out)
        return out




class PropagatorOCbwd(Propagator):
    def __init__(self):
        super().__init__()

    def set_propagator(self, dt, molecule, env):
        self.set_attributes(dt, molecule, env)
        self.clean_propagator()
        self.add_term_to_propagator("eulero1_coeff")
        self.add_term_to_propagator("eulero_energy")
        self.add_term_to_propagator("eulero_field")
        if (self.pcm.env == "sol"):
            self.add_term_to_propagator_bwd("oc_pcm_bwd")


    def propagate_one_step(self, field, wf_fwd):
        for func in self.propagator:
            func(1, field, wf_fwd)

    def propagate_n_step(self, nstep, field, wf_fwd):
        out = list()
        out.append(self.mol.wf.ci)
        for i in range(nstep):
            self.propagate_one_step(field[i], wf_fwd)
            out.append(self.mol.wf.ci)
        out = np.array(out)
        return out



    def add_term_to_propagator_bwd(self, term_name):
        if term_name == "oc_pcm_bwd":
            self.propagator.append(self.bwd_PCM_term)



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






