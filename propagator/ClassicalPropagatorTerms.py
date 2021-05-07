import numpy as np

from propagator.ABCPropagatorTerms import ABCPropagatorTerms
from read_and_set.read import auxiliary_functions as af
from molecule.Molecule import Molecule
from propagator import math_functions as mf

# all posible propagator term. Some use PCM methods specific of FrozenSolventPCM child class, If/when
# DinamicPCM(PCM) will have rabitz implemented, it will have a propagate_bwd_oc term, while the propagation in
# DinamicPCM(PCM) is always done with fortran, so propagate_fortran doesn't exist
# terms are added to the propagator that is used for the wave function propagation
class ClassicalPropagatorTerms(ABCPropagatorTerms):
    def __init__(self):
        super().__init__()
        self.dict_terms = {}


    def init(self):
        self.dict_terms["eulero1_coeff"] = self.eulero1_coeff_term
        self.dict_terms["eulero2_coeff"] = self.eulero2_coeff_term
        self.dict_terms["eulero_energy"] = self.eulero_energy_term
        self.dict_terms["eulero_field"] = self.eulero_field_term
        self.dict_terms["eulero_medium"] = self.eulero_medium_term_fortran
        self.dict_terms["norm"] = self.norm
        self.dict_terms["oc_medium_bwd"] = self.bwd_medium_term_fortran

    # <editor-fold desc="H terms">

    def norm(self, mol, *args):
        mol.wf.norm_ci()

    def eulero1_coeff_term(self, mol, *args):
        mol.wf.ci_prev[0] = np.copy(mol.wf.ci)
        mol.wf.ci = np.copy(mol.wf.ci_prev[0])

    def eulero2_coeff_term(self, mol, *args):
        mol.wf.ci_prev[1] = np.copy(mol.wf.ci_prev[0])
        mol.wf.ci_prev[0] = np.copy(mol.wf.ci)
        mol.wf.ci = np.copy(mol.wf.ci_prev[1])

    def eulero_energy_term(self, mol, order, dt, *args):
        mol.wf.ci += -order * 1j * dt * (mol.par.en_ci * mol.wf.ci_prev[0])

    def eulero_field_term(self, mol, order, dt, field_dt_vector, *args):
        mol.wf.ci += -order * 1j * dt * (-np.dot(np.dot(mol.wf.ci_prev[0], mol.par.muT), field_dt_vector))


    def eulero_medium_term(self, mol, order, dt, field_dt_vector, medium, *args):
        if medium.par.medium != "vac":
            medium.propagate_charges(mol, field_dt_vector)
            mol.wf.ci += -order * 1j * dt \
                          * (np.dot(mol.wf.ci_prev[0],
                                    af.single_summation_tessere(medium.get_q_t(), mol.par.Vijn)))


    def bwd_medium_term(self, mol, order, dt, field_dt_vector, medium, wf_fwd, *args):
        medium.propagate_bwd_oc(wf_fwd, field_dt_vector)
        mol.wf.ci += -order * 1j * dt \
                          * (np.dot(mol.wf.ci_prev[0],
                                    af.single_summation_tessere(medium.get_q_t(), mol.par.Vijn))
                             + np.dot(wf_fwd,
                                      af.single_summation_tessere(
                                            af.double_summation(mol.wf.ci_prev[0],
                                                                np.conj(wf_fwd),
                                                                mol.par.Vijn),
                                            medium.qijn))
                             - np.dot(wf_fwd,
                                      af.single_summation_tessere(
                                            af.double_summation(wf_fwd,
                                                                np.conj(mol.wf.ci_prev[0]),
                                                                mol.par.Vijn),
                                            medium.qijn)))


    def eulero_medium_term_fortran(self, mol, order, dt, field_dt_vector, medium, *args):
        if medium.par.medium != "vac":
            medium.propagate_charges_fortran(mol, field_dt_vector)
            q_t = medium.get_q_t()
            mol.wf.ci += -order * 1j * dt \
                      *(mf.eulero_pcm(np.asfortranarray(mol.wf.ci_prev[0], dtype=np.complex128),
                                      np.asfortranarray(q_t, dtype=np.complex128),
                                      np.asfortranarray(mol.par.Vijn_fortran_flip, dtype=np.complex128)))



    def bwd_medium_term_fortran(self, mol, order, dt, field_dt_vector, medium, wf_fwd, *args):
        medium.propagate_bwd_oc_fortran(wf_fwd, field_dt_vector)
        q_t = medium.get_q_t()
        mol.wf.ci += -order * 1j * dt \
                          * (mf.bwd_pcm(np.asfortranarray(mol.wf.ci_prev[0], dtype=np.complex128),
                                        np.asfortranarray(wf_fwd, dtype=np.complex128),
                                        np.asfortranarray(q_t, dtype=np.complex128),
                                        np.asfortranarray(mol.par.Vijn_fortran_flip, dtype=np.complex128),
                                        np.asfortranarray(medium.qijn_fortran_flip, dtype=np.complex128)))



