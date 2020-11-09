# -*- coding: utf-8 -*-

import dictionaries.OCDictionaries as ocdict
from molecule.Molecule import Molecule
from ObjectiveFunctionTerms import ObjectiveFunctionTerms


class ObjectiveFunction():
    def __init__(self):
        self.objective_function = []
        self.function_terms = ObjectiveFunctionTerms()
        self.observable = None
        
        
    def clean_objective_function(self):
        self.objective_function = []
        
        
    def add_term_to_obj(self, term_name):
        self.objective_function.append(self.function_terms.dict_terms[term_name])
        
    
    def compute_objective_function(self, *args):
        pass


