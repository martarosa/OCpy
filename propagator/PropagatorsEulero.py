import numpy as np
from propagator.ABCPropagator import ABCPropagator
from propagator.PropagatorTerms import PropagatorTerms

class PropagatorEulero1Order(ABCPropagator):
    def __init__(self):
        super().__init__()
        self.propagator_terms = PropagatorTerms()
        self.propagator = []

    def set_propagator(self, molecule, env):
        self.init_propagaror_terms(molecule, env)
        self.clean_propagator()
        self.add_term_to_propagator("eulero1_coeff")
        self.add_term_to_propagator("eulero_energy")
        self.add_term_to_propagator("eulero_field")
        if self.propagator_terms.pcm != None:
            if self.propagator_terms.pcm.par.env == "sol":
                self.add_term_to_propagator("eulero_pcm")
        self.add_term_to_propagator("norm")

    def propagate_one_step(self, i, dt, field):
        for func in self.propagator:
            func(i, 1, dt, field)

    def propagate_n_step(self, dt, nstep, field):
        out = list()
        out.append(self.propagator_terms.mol.wf.ci)
        for i in range(nstep):
            self.propagate_one_step(i, dt, field[i])
            out.append(self.propagator_terms.mol.wf.ci)
        out = np.array(out)
        return out



class PropagatorEulero2Order(ABCPropagator):
    def __init__(self):
        super().__init__()
        self.propagator_terms = PropagatorTerms()
        self.propagator = []

    def set_propagator(self, molecule, env):
        self.init_propagaror_terms(molecule, env)
        self.clean_propagator()
        self.add_term_to_propagator("eulero2_coeff")
        self.add_term_to_propagator("eulero_energy")
        self.add_term_to_propagator("eulero_field")
        if self.propagator_terms.pcm != None:
            if self.propagator_terms.pcm.par.env == "sol":
                self.add_term_to_propagator("eulero_pcm")
        self.add_term_to_propagator("norm")


    def propagate_one_step(self, i, dt, field, order=2):
        for func in self.propagator:
            func(i, order, dt, field)

    def propagate_n_step(self, dt, nstep, field):
        out = list()
        out.append(self.propagator_terms.mol.wf.ci)
        for i in range(nstep):
            if i != 0:
                self.propagate_one_step(i, dt, field[i])
            else:
                self.propagate_one_step(i, dt, field[i], order=1)
            out.append(self.propagator_terms.mol.wf.ci)
        out = np.array(out)
        return out