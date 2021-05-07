from abc import ABCMeta, abstractmethod

from propagator.ABCPropagatorTerms import ABCPropagatorTerms
from molecule.Molecule import Molecule
from medium.ABCMedium import ABCMedium
from dictionaries import Dictionaries as dict


#Propagator contains Molecule() and Medium() objects, which are stored only here.
#In Medium() there are the methods to propagate the charges, that must be called from there
#Here we decided to use "delegates" to perform he propagation, here explained:
#In propagator terms there is a different method for each term of the propagator (Hamiltonian). The methods receive
#the information needed, particularly the time dependent ones (c_t, q_t, field_t) and gives back c_t+1 for each
#single propagation term.
#The needed terms are added in the self.propagator list with the set_propagator() method (which is different for each
#child of ABCPropagator
#propagate_one_step and propagate_n_step run though all the element of the list.
#Each of them sum to the c_t+1 which results from the application of the previous term to c_t, the c_t+1 resulting from
#its own calculation


class ABCPropagator(metaclass=ABCMeta):
    def __init__(self):
        self.mol = Molecule()
        self.medium = None #ABCMedium()

        self.propagator_terms = ABCPropagatorTerms()
        self.propagator = []


    def init(self, molecule, medium):
        self.mol = molecule
        self.medium = medium
        self.propagator_terms = dict.PropagatorTermsDict["classical"]()
        self.propagator_terms.init()



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








