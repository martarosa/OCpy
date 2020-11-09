# -*- coding: utf-8 -*-

import numpy as np
from read_and_set.read import auxiliary_functions as af

class ObjectiveFunctionTerms():
    def __init__(self):
        self.dict_terms = {}
        self.propagation_type = None
        
        
    def init(self):
        if self.propagation_type == "trotter_suzuki":
            self.dict_terms["simple_observable"] == self.expectation_value_observable_quantum
        else: 
            self.dict_terms["simple_observable"] = self.expectation_value_observable
        self.dict_terms["field_fluency_penalization"] = self.field_fluency_penalization
    

### siamo sicuri che come argomento ci vada oc?
    
    def expectation_value_observable(self, oc, observable):
        oc.par.J = af.compute_expectation_value(oc.mol.wf.ci, observable)
        
        
    def expectation_value_observable_quantum(self):
        pass
        
        
    def field_fluency_penalization(self, oc):
        oc.par.J -= oc.alpha_field_J_integral() 
        
