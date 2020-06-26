from propagator.PropagatorTerms import PropagatorTerms
from abc import ABCMeta, abstractmethod


# in propagator_terms there are all possible propagation terms (probably we could delete this attribute and
# temporary inlitialized in each child class.
# Then in each child class the specific terms arre added to the propagator delegate, and the propagation is done
#cycling through all the terms in propagate_one_step (for funct in propagator: func(..))




class ABCPropagator(metaclass=ABCMeta):
    def __init__(self):
        self.propagator_terms = PropagatorTerms()
        self.propagator = []



    def init_propagator_terms(self, molecule, medium):
        self.propagator_terms.set_attributes(molecule, medium)
        self.propagator_terms.init_terms_dictionary()

    def add_term_to_propagator(self, term_name):
            self.propagator.append(self.propagator_terms.dict_terms[term_name])

    def clean_propagator(self):
        self.propagator = []

    @abstractmethod
    def set_propagator(self, molecule, medium):
        pass

    @abstractmethod
    #takes vector and matricise (field) relative to the specific time of the propagation
    #and gives back a n_states vector of the propagated wavefunction
    def propagate_one_step(self, *args):
        pass

    @abstractmethod
    #takest objects of tipe Func_t with a vector with values of the time axes and a matrix of the values in each time step
    #e.g. field is fx, fy,fz and returns an object Func_t for the wavefuncion, with the time axis and a matrix of the
    # wf in all the states
    def propagate_n_step(self, *args):
        pass











