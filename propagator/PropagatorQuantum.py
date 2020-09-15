#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  4 12:04:17 2020

@author: davide
"""


from qiskit import Aer, execute
from medium.ABCMedium import ABCMedium
from molecule.Molecule import Molecule
from propagator.ABCPropagator import ABCPropagator
from propagator.QuantumPropagatorTermsIBM import QuantumPropagatorTermsIBM



class PropagatorQuantum(ABCPropagator):
    def __init__(self):
        super().__init__()
        self.mol = Molecule()
        self.medium = ABCMedium()
        self.propagator_terms = QuantumPropagatorTermsIBM()
        self.propagator = []


    def set_propagator(self, molecule, medium):
        self.mol = molecule
        self.medium = medium
        self.clean_propagator()
        self.propagator_terms.init(molecule)
        self.add_term_to_propagator("quantum_evolution")

    def propagate_one_step(self, dt, field_dt_vector):
        for func in self.propagator:
            func(dt, field_dt_vector)

    def propagate_n_step(self, discrete_time_par, field):
        for i in range(discrete_time_par.nstep-1):
            field_dt_vector = field.f_xyz[i] + field.f_xyz[i+1]
            self.propagate_one_step(discrete_time_par.dt, field_dt_vector)
        result = execute(self.propagator_terms.qcircuit, self.propagator_terms.provider, shots=1).result()
        statevector = result.get_statevector(self.propagator_terms.qcircuit) #da sistemare
        return statevector