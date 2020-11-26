# -*- coding: utf-8 -*-

from ObjectiveFunction.ObjectiveFunctionTerms import ObjectiveFunctionTerms
from ObjectiveFunction.ABCObjectiveFunction import ABCObjectiveFunction


class ObjectiveFunctionOptical(ABCObjectiveFunction):
    def __init__(self):
        super().__init__()
        
    
    def set_objective_function(self):
        self.function_terms.init()
        self.add_term_to_obj("target_state_population")
        self.add_term_to_obj("field_fluency_penalization")
    

    def compute_objective_function(self, oc, wavefunction, target_state):
        for func in self.objective_function:
            func(oc, wavefunction, target_state)
    
    
class ObjectiveFunctionGroundState(ABCObjectiveFunction):
    def __init__(self):
        super().__init__()
        
    
    def set_objective_function(self):
        self.function_terms.init()
        self.add_term_to_obj("simple_observable")
    

    def compute_objective_function(self, oc, wavefunction, observable):  
        for func in self.objective_function:
            func(oc, wavefunction, observable)
        
    
    