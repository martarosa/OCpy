#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  6 12:04:51 2020

@author: davide
"""

import numpy as np

from read_and_set.read import auxiliary_functions as af
from molecule.Molecule import Molecule
#from propagator import math_functions as mf

from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister
from propagator.PropagatorTerms import PropagatorTerms


from read_and_set.read import NamelistTools as nmtool
import sys


class QuantumPropagatorTerms(PropagatorTerms):
        def __init__(self):
            self.mol = Molecule()
            self.pcm = None
            self.dict_terms = {}
            self.is_quantum = None
            self.qbits = None
            self.cbits = None
            self.qcircuit = None
            
            
            
        def set_qbits(self):
            self.qbits = QuantumRegister(self.mol.wf.n_ci, 'q')
            
        def set_cbits(self):
            self.cbits = ClassicalRegister(self.mol.wf.n_ci, 'c')
            
        def set_qcircuit(self):
            self.qcircuit = QuantumCircuit(self.qbits,self.cbits)
            
                
        def set_attributes(self, molecule, pcm):
            super().set_attributes(molecule, pcm)
            self.qbits = QuantumRegister(self.mol.wf.n_ci, 'q')
            self.cbits = ClassicalRegister(self.mol.wf.n_ci, 'c')
            self.qcircuit = QuantumCircuit(self.qbits,self.cbits)

        def init_terms_dictionary(self):
            self.dict_terms["quantum_evo"] = self.quantum_evo_circuit



        def quantum_evo_circuit(self, i, dt, field_dt_vector):
            if i == 0:
                self.qcircuit.x(self.qbits[0])
            for k in range(self.mol.wf.n_ci):
                self.qcircuit.u1((self.mol.par.en_ci[k] - np.dot(self.mol.par.muT,(field_dt_vector)/2)[k][k])*dt, self.qbits[k])
                self.qcircuit.h(self.qbits[k])
            self.bounded_occupation_gatex(dt, field_dt_vector)
            for k in range(self.mol.wf.n_ci):
                self.qcircuit.rx(np.pi/2, self.qbits[k])
            self.bounded_occupation_gatey(dt, field_dt_vector)
            
        def bounded_occupation_gatex(self, dt, field_dt_vector):
            for i in range(self.mol.wf.n_ci-1):
                for j in range(1,self.mol.wf.n_ci):
                    if j > i:
                        self.qcircuit.cx(self.qbits[j], self.qbits[i])
                        self.qcircuit.rz(np.dot(self.mol.par.muT,((field_dt_vector)/2))[j][i]*dt, self.qbits[i])
                        self.qcircuit.cx(self.qbits[j], self.qbits[i])
                self.qcircuit.h(self.qbits[i])
            self.qcircuit.h(self.qbits[j])
            
        def bounded_occupation_gatey(self, dt, field_dt_vector):
            for i in range(self.mol.wf.n_ci-1):
                for j in range(1,self.mol.wf.n_ci):
                    if j > i:
                        self.qcircuit.cx(self.qbits[j], self.qbits[i])
                        self.qcircuit.rz(np.dot(self.mol.par.muT,((field_dt_vector)/2))[j][i]*dt, self.qbits[i])
                        self.qcircuit.cx(self.qbits[j], self.qbits[i])
                self.qcircuit.rx(-np.pi/2, self.qbits[i])
            self.qcircuit.rx(-np.pi/2, self.qbits[j])
            