import numpy as np

from propagator import PropagatorsEulero as prop
from field.Field import Field
from read import auxiliary_functions as af

class Chromosome():
    def __init__(self):
        self.size_amplitude_list = None
        self.amplitudes = []  # n_ampitudes
        self.J = None

        self.field = Field()
        self.prop_psi = prop.PropagatorEulero2Order()


    def init_chromosome(self, dt, molecule, field, pcm):
        self.field = field
        self.size_amplitude_list = int(self.field.parameters['fi'].size*2)
        self.prop_psi.set_propagator(dt, molecule, pcm)
        self.field_to_amplitudes_list()


    def field_to_amplitudes_list(self):
        self.amplitudes = np.concatenate([self.field.parameters['fi'].ravel(),

                                          self.field.parameters['fi_cos'].ravel()]).tolist()

    def amplitudes_to_field(self):
        self.field.parameters['fi'] = np.asarray(self.amplitudes[:int(self.size_amplitude_list/2)]).reshape((-1,3))
        self.field.parameters['fi_cos'] = np.asarray(self.amplitudes[:int(self.size_amplitude_list / 2)]).reshape((-1, 3))
        self.field.chose_field('sum')



    def calc_J(self, target_state, alpha_t, dt):
        self.J = np.real(af.projector_mean_value(
            self.prop_psi.propagator_terms.mol.wf.ci,
            target_state)
                         - af.alpha_field_J_integral(
            self.field.field,
            alpha_t, dt))