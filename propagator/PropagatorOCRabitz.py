import numpy as np

from propagator.ABCPropagator import ABCPropagator
from propagator.PropagatorTerms import PropagatorTerms


#import math_functions as mf

class PropagatorOCfwd(ABCPropagator):
    def __init__(self):
        super().__init__()
        self.propagator_terms = PropagatorTerms()
        self.propagator = []

    def set_propagator(self, dt, molecule, env):
        self.init_propagaror_terms(dt, molecule, env)
        self.clean_propagator()
        self.add_term_to_propagator("eulero1_coeff")
        self.add_term_to_propagator("eulero_energy")
        self.add_term_to_propagator("eulero_field")
        if self.propagator_terms.pcm.env == "sol":
            self.add_term_to_propagator("eulero_pcm")


    def propagate_one_step(self, i, field):
        for func in self.propagator:
            func(i, 1, field)


    def propagate_n_step(self, nstep, field):
        out = list()
        out.append(self.propagator_terms.mol.wf.ci)
        for i in range(nstep):
            self.propagate_one_step(i, field[i])
            out.append(self.propagator_terms.mol.wf.ci)
        out = np.array(out)
        return out


class PropagatorOCbwd(ABCPropagator):
    def __init__(self):
        super().__init__()
        self.propagator_terms = PropagatorTerms()
        self.propagator = []

    def set_propagator(self, dt, molecule, env):
        self.init_propagaror_terms(dt, molecule, env)
        self.clean_propagator()
        self.add_term_to_propagator("eulero1_coeff")
        self.add_term_to_propagator("eulero_energy")
        self.add_term_to_propagator("eulero_field")
        if (self.propagator_terms.pcm.env == "sol"):
            self.add_term_to_propagator("oc_pcm_bwd")


    def propagate_one_step(self, i, field, wf_fwd):
        for func in self.propagator:
            func(i, 1, field, wf_fwd)

    def propagate_n_step(self, nstep, field, wf_fwd):
        out = list()
        out.append(self.propagator_terms.mol.wf.ci)
        for i in range(nstep):
            self.propagate_one_step(i, field[i], wf_fwd)
            out.append(self.propagator_terms.mol.wf.ci)
        out = np.array(out)
        return out










