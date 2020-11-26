# -*- coding: utf-8 -*-

from abc import ABCMeta, abstractmethod
from ObjectiveFunction.ObjectiveFunctionTerms import ObjectiveFunctionTerms


class ABCObjectiveFunction(metaclass=ABCMeta):
    def __init__(self):
        self.objective_function = []
        self.function_terms = ObjectiveFunctionTerms()
        self.observable = None
        self.J = None
        
        
    def clean_objective_function(self):
        self.objective_function = []
        
        
    def add_term_to_obj(self, term_name):
        self.objective_function.append(self.function_terms.dict_terms[term_name])
        
    @abstractmethod
    def set_objective_function():
        pass
        
    @abstractmethod
    def compute_objective_function(self, *args):
        pass


