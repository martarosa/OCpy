from abc import ABCMeta, abstractmethod

import dictionaries.PropagatorDictionaries
from propagator.ABCPropagatorTerms import ABCPropagatorTerms
from molecule.Molecule import Molecule
from medium.ABCMedium import ABCMedium
from dictionaries import SaveDictionaries as dict

class ABCPropagator(metaclass=ABCMeta):
    def __init__(self):
        self.mol = Molecule()
        self.medium = None

        self.propagator_terms = ABCPropagatorTerms()
        self.propagator = []


    def init(self, molecule, medium, prop_conf, propagator = None):
        self.mol = molecule
        self.medium = medium
        self.propagator_terms = dictionaries.PropagatorDictionaries.PropagatorTermsDict[propagator]()
        self.propagator_terms.init()


    def add_term_to_propagator(self, term_name):
        self.propagator.append(self.propagator_terms.dict_terms[term_name])

    def clean_propagator(self):
        self.propagator = []

    @abstractmethod
    def set_propagator(self, molecule, medium, prop_conf):
        pass

    @abstractmethod
    #takes vector and matricise (field) relative to the specific time of the propagation
    #and gives back a n_states vector of the propagated wavefunction
    def propagate_one_step(self, *args):
        pass

    @abstractmethod
    #takes objects of tipe Func_t with a vector with values of the time axes and a matrix of the values in each time step
    #e.g. field is fx, fy,fz and returns an object Func_t for the wavefuncion, with the time axis and a matrix of the
    # wf in all the states
    def propagate_n_step(self, *args):
        pass








