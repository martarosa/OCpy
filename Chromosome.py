import numpy as np
from copy import deepcopy

from propagator import PropagatorsEulero as prop
from field.Field import Field
from read_and_set.read import auxiliary_functions as af

class Chromosome():
    def __init__(self):
        self.n_ampl = None
        self.amplitudes = []  # n_ampitudes
        self.J = None
        self.field = Field()
        self.prop_psi = prop.PropagatorEulero2Order()


    def init_chromosome(self, dt, molecule, field, pcm):
        self.field = field
        self.n_ampl = int(self.field.parameters['fi'].size)
        self.prop_psi.set_propagator(dt, molecule, pcm)
        self.amplitudes = self.field.parameters['fi'].ravel().tolist()


 #       self.field_to_amplitudes()


#    def field_to_amplitudes(self):
#        self.amplitudes = self.field.parameters['fi'].ravel().tolist()


    def amplitudes_to_field(self):
        self.field.parameters['fi'] = np.asarray(self.amplitudes).reshape((-1,3))
        self.field.chose_field('sum')


    def calc_J(self, target_state, alpha_t, dt):
        Jprev = deepcopy(self.J)
        self.J = np.real(af.projector_mean_value(
                                       self.prop_psi.propagator_terms.mol.wf.ci,
                                       target_state)
                         -self.alpha_field_J_integral(
                                       self.field.fieldaaaaaaaaaaaaaaaaaaa,
                                       alpha_t,
                                       dt))
        if self.J > Jprev:
            return 1
        else:
            return 0