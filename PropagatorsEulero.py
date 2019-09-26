import numpy as np
from Propagator import Propagator


class PropagatorEulero1Order(Propagator):
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
        self.add_term_to_propagator("norm")

    def propagate_one_step(self, field):
        for func in self.propagator:
            func(1,field)

    def propagate_n_step(self, nstep, field):
        out = list()
        out.append(self.mol.wf.ci)
        for i in range(nstep):
            self.propagate_one_step(field[i])
            out.append(self.mol.wf.ci)
        out = np.array(out)
        return out



class PropagatorEulero2Order(Propagator):
    def __init__(self):
        super().__init__()

    def set_propagator(self, dt, molecule, env):
        self.set_attributes(dt, molecule, env)
        self.clean_propagator()
        self.add_term_to_propagator("eulero2_coeff")
        self.add_term_to_propagator("eulero_energy")
        self.add_term_to_propagator("eulero_field")
        if self.pcm.env == "sol":
            self.add_term_to_propagator("eulero_pcm")
        self.add_term_to_propagator("norm")


    def propagate_one_step(self, field, order=2):
        for func in self.propagator:
            func(order, field)

    def propagate_n_step(self, nstep, field):
        out = list()
        out.append(self.mol.wf.ci)
        for i in range(nstep):
            if i != 0:
                self.propagate_one_step(field[i])
            else:
                self.propagate_one_step(field[i], order=1)
            out.append(self.mol.wf.ci)
        out = np.array(out)
        return out