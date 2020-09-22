# -*- coding: utf-8 -*-
"""
Created on Tue Sep 22 11:22:48 2020

@author: castd
"""

import numpy as np


from propagator.ABCPropagator import ABCPropagator
from propagator.QuantumPropagatorTerms import QuantumPropagatorTerms


class PropagatorQuantum(ABCPropagator):
    def __init__(self):
        super().__init__()    ### perchè scrivere super().__init__() se dopo richiamo tutti i metodi?
        self.propagator_terms = QuantumPropagatorTerms()
        self.propagator = []

    def set_propagator(self, molecule, medium):
        self.init(molecule, medium)
        self.clean_propagator()
        self.add_term_to_propagator("quantum_evo")
        
    def propagate_one_step(self, i, dt, field_dt_vector):
        for func in self.propagator:
            func(i, dt, field_dt_vector)
    
    def propagate_n_step(self, discrete_time_par, field):
        for i in range(discrete_time_par.nstep-1):
            field_dt_vector = field.f_xyz[i] + field.f_xyz[i+1]
            self.propagate_one_step(i, discrete_time_par.dt, field_dt_vector)
        result = self.propagator_terms.execute()
        return result   