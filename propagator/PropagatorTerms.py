import numpy as np

from read_and_set.read import auxiliary_functions as af
from molecule.Molecule import Molecule


# all posible propagator term. Some use PCM methods specific of FrozenSolventPCM child class, If/when
# DinamicPCM(PCM) will have rabitz implemented, it will have a propagate_bwd_oc term, while the propagation in
# DinamicPCM(PCM) is always done with fortran, so propagate_fortran doesn't exist
# terms are added to the propagator that is used for the wave function propagation
class PropagatorTerms():
    def __init__(self):
        self.mol = Molecule()
        self.pcm = None #ABCPCM()
        self.dict_terms = {}


    def set_attributes(self, molecule, pcm):
        self.mol = molecule
        self.pcm = pcm

    def init_terms_dictionary(self):
        self.dict_terms["eulero1_coeff"] = self.eulero1_coeff_term
        self.dict_terms["eulero2_coeff"] = self.eulero2_coeff_term
        self.dict_terms["eulero_energy"] = self.eulero_energy_term
        self.dict_terms["eulero_field"] = self.eulero_field_term
        self.dict_terms["eulero_pcm"] = self.eulero_PCM_term
        self.dict_terms["norm"] = self.norm
        self.dict_terms["oc_pcm_bwd"] = self.bwd_PCM_term

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

    def eulero_energy_term(self, i, order, dt, *args):
        self.mol.wf.ci += -order * 1j * dt * (self.mol.par.en_ci * self.mol.wf.ci_prev[0])

    def eulero_field_term(self, i, order, dt, field_dt_vector, *args):
        self.mol.wf.ci += -order * 1j * dt * (-np.dot(np.dot(self.mol.wf.ci_prev[0], self.mol.par.muT), field_dt_vector))

    def eulero_PCM_term(self, i, order, dt, field_dt_vector, *args):
        self.pcm.propagate(i, self.mol, field_dt_vector)
        self.mol.wf.ci += -order * 1j * dt \
                          * (np.dot(self.mol.wf.ci_prev[0],
                                    af.single_summation_tessere(self.pcm.get_q_t() - self.pcm.q00n, self.mol.par.Vijn)))

    def bwd_PCM_term(self, i, order, dt, field_dt_vector, wf_fwd, *args):
        self.pcm.propagate_bwd_oc(i, self.mol, field_dt_vector, wf_fwd)
        self.mol.wf.ci += -order * 1j * dt \
                          * (np.dot(self.mol.wf.ci_prev[0],
                                    af.single_summation_tessere(self.pcm.get_q_t() - self.pcm.q00n,
                                                                self.mol.par.Vijn))
                             + np.dot(wf_fwd,
                                      af.single_summation_tessere(
                                            af.double_summation(self.mol.wf.ci_prev[0],
                                                                np.conj(wf_fwd),
                                                                self.mol.par.Vijn),
                                            self.pcm.qijn))
                             - np.dot(wf_fwd,
                                      af.single_summation_tessere(
                                            af.double_summation(wf_fwd,
                                                                np.conj(self.mol.wf.ci_prev[0]),
                                                                self.mol.par.Vijn),
                                            self.pcm.qijn)))







