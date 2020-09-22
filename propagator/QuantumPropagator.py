# -*- coding: utf-8 -*-
"""
Created on Tue Sep 22 11:22:48 2020

@author: castd
"""

import numpy as np
from qiskit import QuantumRegister, ClassicalRegister, QuantumCircuit
from qiskit import Aer, execute

from propagator.ABCPropagator import ABCPropagator
from propagator.QuantumPropagatorTerms import QuantumPropagatorTerms


class PropagatorQuantum(ABCPropagator):
    def __init__(self):
        super().__init__()    ### perchè scrivere super().__init__() se dopo richiamo tutti i metodi?
        self.propagator_terms = QuantumPropagatorTerms()
        self.propagator = []
        self.provider = None
        self.qbits = None
        self.cbits = None
        self.qcircuit = None


    def set_qprocessor(self):
        self.propagator_terms.set_qbits()
        self.propagator_terms.set_cbits()
        self.propagator_terms.set_qcircuit()
        
    def set_propagator(self, molecule, env):
        self.init_propagaror_terms(molecule, env)
        self.clean_propagator()
        self.set_qprocessor()
        self.add_term_to_propagator("quantum_evo")
        
    def propagate_one_step(self, i, dt, field_dt_vector):
        for func in self.propagator:
            func(i, dt, field_dt_vector)
    
    def propagate_n_step(self, discrete_time_par, field):
        for i in range(discrete_time_par.nstep-1):
            field_dt_vector = field.f_xyz[i] + field.f_xyz[i+1]
            self.propagate_one_step(i, discrete_time_par.dt, field_dt_vector)
        result = execute(self.propagator_terms.qcircuit, self.provider, shots=1).result()
        statevector = result.get_statevector(self.propagator_terms.qcircuit)
        return statevector   #### non farti restituire statevector crea un attributo da aggiornare o forse lascia cosi, propagate_nstep classico restituisce tutta la traiettoria