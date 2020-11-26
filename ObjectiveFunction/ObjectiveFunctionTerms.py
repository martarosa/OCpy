# -*- coding: utf-8 -*-

import numpy as np
from read_and_set.read import auxiliary_functions as af

class ObjectiveFunctionTerms():
    def __init__(self):
        self.dict_terms = {}
        self.propagation_type = None
        
        
    def init(self):
        if self.propagation_type == "trotter_suzuki":
            self.dict_terms["simple_observable"] = self.expectation_value_observable_quantum
            self.dict_terms["target_state_population"] = self.target_state_population_quantum
        else: 
            self.dict_terms["simple_observable"] = self.expectation_value_observable
            self.dict_terms["target_state_population"] = self.target_state_population
        self.dict_terms["field_fluency_penalization"] = self.field_fluency_penalization
    

### siamo sicuri che come argomento ci vada oc?
    
    def expectation_value_observable(self, oc, wavefunction, observable):
        oc.obj_fun.J = af.compute_expectation_value(wavefunction, observable)
        
        
    def expectation_value_observable_quantum(self):
        pass
        
        
    def field_fluency_penalization(self, oc):
        oc.obj_fun.J -= oc.alpha_field_J_integral()
        
        
    def target_state_population(self, oc, wavefunction, target_state):
        oc.obj_fun.J = af.projector_mean_value(wavefunction, target_state)
        
        
    def target_state_population_quantum(self, oc):
        pass
        
