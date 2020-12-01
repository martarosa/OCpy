# -*- coding: utf-8 -*-
"""
Created on Tue Sep 22 11:22:48 2020

@author: castd
"""

import numpy as np


from propagator.ABCPropagator import ABCPropagator
from propagator.QuantumPropagatorTerms import QuantumPropagatorTerms
import dictionaries.PropagatorTermsDictionaries as ptdict

from molecule.Molecule import Molecule


class PropagatorQuantum(ABCPropagator):
    def __init__(self):
        self.mol = Molecule()
        self.medium = None  
        self.propagator_terms = None
        self.propagator = []
        self.final_state_qc = None
        self.counts_dictionary = None

    def set_propagator(self, molecule, medium, prop_conf):
        self.init(molecule, medium, prop_conf)
        self.clean_propagator(molecule)
        self.add_term_to_propagator("quantum_evo")
        
        
    def clean_propagator(self, mol):
        self.propagator = []
        self.propagator_terms.set_qprocessor(mol)
        
    def init(self, molecule, medium, propagator):
        self.mol = molecule
        self.medium = medium
        print(propagator.quantum_prop_keyword)
        self.propagator_terms = ptdict.PropagatorTermsDict[propagator.quantum_prop_keyword]()   
        self.propagator_terms.init(propagator)
        
    def propagate_one_step(self, dt, field_dt_vector):
        for func in self.propagator:
            func(self.mol, dt, field_dt_vector)
    
    def propagate_n_step(self, discrete_time_par, field):
        for i in range(discrete_time_par.nstep-1):
            if i == 0:
                self.propagator_terms.prepare_GS_linear_mapping()
            field_dt_vector = field.f_xyz[i] + field.f_xyz[i+1]
            self.propagate_one_step(discrete_time_par.dt, field_dt_vector)
            if i == discrete_time_par.nstep-1 and self.propagator_terms.IBMInterface.provider != "statevector_simulator":
                self.propagator_terms.computational_basis_measurement()
        result = self.propagator_terms.execute()
        if type(result[0]) == float:
            self.counts_dictionary = result
        else:
            self.final_state_qc = result
    
        
    