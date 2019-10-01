import numpy as np
from propagator.Propagator import Propagator
from propagator.PropagatorTerms import PropagatorTerms

class PropagatorEulero1Order(Propagator):
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
        self.add_term_to_propagator("norm")

    def propagate_one_step(self, field):
        for func in self.propagator:
            func(1,field)

    def propagate_n_step(self, nstep, field):
        out = list()
        out.append(self.propagator_terms.mol.wf.ci)
        for i in range(nstep):
            self.propagate_one_step(field[i])
            out.append(self.propagator_terms.mol.wf.ci)
        out = np.array(out)
        return out



class PropagatorEulero2Order(Propagator):
    def __init__(self):
        super().__init__()
        self.propagator_terms = PropagatorTerms()
        self.propagator = []

    def set_propagator(self, dt, molecule, env):
        self.init_propagaror_terms(dt, molecule, env)
        self.clean_propagator()
        self.add_term_to_propagator("eulero2_coeff")
        self.add_term_to_propagator("eulero_energy")
        self.add_term_to_propagator("eulero_field")
        if self.propagator_terms.pcm.env == "sol":
            self.add_term_to_propagator("eulero_pcm")
        self.add_term_to_propagator("norm")


    def propagate_one_step(self, field, order=2):
        for func in self.propagator:
            func(order, field)

    def propagate_n_step(self, nstep, field):
        out = list()
        out.append(self.propagator_terms.mol.wf.ci)
        for i in range(nstep):
            if i != 0:
                self.propagate_one_step(field[i])
            else:
                self.propagate_one_step(field[i], order=1)
            out.append(self.propagator_terms.mol.wf.ci)
        out = np.array(out)
        return out