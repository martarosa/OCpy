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

from qiskit import Aer, execute
from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister

from propagator.ABCPropagatorTerms import ABCPropagatorTerms


from read_and_set.read import NamelistTools as nmtool
import sys


class QuantumPropagatorTermsIBM(ABCPropagatorTerms):
        def __init__(self):
            super().__init__()
            self.dict_terms = {}

            self.qbits = None
            self.cbits = None
            self.qcircuit = None
            self.provider = None

        def init(self, molecule):
            self.qbits = QuantumRegister(molecule.wf.n_ci, 'q')
            self.cbits = ClassicalRegister(molecule.wf.n_ci, 'c')
            self.qcircuit = QuantumCircuit(self.qbits,self.cbits)
            self.qcircuit.x(self.qbits[0])
            self.init_terms()


        def init_terms(self):
            self.dict_terms["quantum_evolution"] = self.quantum_evolution_circuit


        def quantum_evolution_circuit(self, dt, mol, field_dt_sum_interval_vector):
            for k in range(mol.wf.n_ci):
                self.qcircuit.u1((mol.par.en_ci[k] - np.dot(mol.par.muT, (field_dt_sum_interval_vector) / 2)[k][k]) * dt, self.qbits[k])
                self.qcircuit.h(self.qbits[k])
            self.bounded_occupation_gatex(dt, mol, field_dt_sum_interval_vector)
            for k in range(mol.wf.n_ci):
                self.qcircuit.rx(np.pi/2, self.qbits[k])
            self.bounded_occupation_gatey(dt, mol, field_dt_sum_interval_vector)

        def bounded_occupation_gatex(self, dt, mol, field_dt_sum_interval_vector):
            for i in range(mol.wf.n_ci-1):
                for j in range(1,mol.wf.n_ci):
                    if j > i:
                        self.qcircuit.cx(self.qbits[j], self.qbits[i])
                        self.qcircuit.rz(np.dot(mol.par.muT, ((field_dt_sum_interval_vector) / 2))[j][i] * dt, self.qbits[i])
                        self.qcircuit.cx(self.qbits[j], self.qbits[i])
                self.qcircuit.h(self.qbits[i])
            self.qcircuit.h(self.qbits[mol.wf.n_ci])

        def bounded_occupation_gatey(self, dt, mol, field_dt_sum_interval_vector):
            for i in range(mol.wf.n_ci-1):
                for j in range(1,mol.wf.n_ci):
                    if j > i:
                        self.qcircuit.cx(self.qbits[j], self.qbits[i])
                        self.qcircuit.rz(np.dot(mol.par.muT, ((field_dt_sum_interval_vector) / 2))[j][i] * dt, self.qbits[i])
                        self.qcircuit.cx(self.qbits[j], self.qbits[i])
                self.qcircuit.rx(-np.pi/2, self.qbits[i])
            self.qcircuit.rx(-np.pi/2, self.qbits[mol.wf.n_ci])
