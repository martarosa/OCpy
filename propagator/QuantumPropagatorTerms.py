# -*- coding: utf-8 -*-
"""
Created on Tue Sep 22 11:03:33 2020

@author: castd
"""

import numpy as np

from read_and_set.read import auxiliary_functions as af
from molecule.Molecule import Molecule
#from propagator import math_functions as mf

from qiskit import QuantumCircuit , QuantumRegister , ClassicalRegister
from qiskit import Aer , execute , IBMQ
from propagator.ABCPropagatorTerms import ABCPropagatorTerms
from qiskit.providers.aer.noise import NoiseModel

import sys


class QuantumPropagatorTerms(ABCPropagatorTerms):
        def __init__(self):
            self.dict_terms = {}
            self.qbits = None
            self.cbits = None
            self.qcircuit = None
            self.IBMInterface = None

            
        def set_qbits(self):
            self.qbits = QuantumRegister(self.mol.wf.n_ci, 'q')
            
        def set_cbits(self):
            self.cbits = ClassicalRegister(self.mol.wf.n_ci, 'c')
            
        def set_qcircuit(self):
            self.qcircuit = QuantumCircuit(self.qbits,self.cbits)
            
            
        def set_qprocessor(self):
            self.propagator_terms.set_qbits()
            self.propagator_terms.set_cbits()
            self.propagator_terms.set_qcircuit()

        def init(self):
            self.q_processor()
            self.dict_terms["quantum_evo"] = self.quantum_evo_circuit
            self.dict_terms["expectation_value"] = self.expectation_value_hermitian_operator
            self.dict_terms["measure_computational_basis"] = self.computational_basis_measurement
                    
            
            
##### Method to build the circuit equivalent to simulation of a single_particle hamiltonian, in the field-free diagonal basis, under a TD external perturbation ##### 


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
            
            
###### Measurement_methods #######
            
        def computational_basis_measurement(self):
            self.qcircuit.measure(self.qbits, self.cbits)
            
        def expectation_value_hermitian_operator(self):
            """ è verosimile che questo metodo avrà dentro l'inizializzazione di vari circuiti con questo schema:
                CIRCUITO DI EVOLUZIONE + STRINGHE DI PAULI VARIE + MISURA BASE COMPUTAZIONALE
                fino ad esaurimento delle stringhe di pauli.
                Visto che execute è separato deve restituire una lista di circuiti da eseguire
                Per ottenere il valore di aspettazione poi avremo una cosa del tipo:
                for counts in data:
                    expectation_value += counts*matrix_element_giusto """
            pass
                     
###### different kinds of execution will be needed #######
            
        def configure_backend(self, simulator_options):  ## simulator_options è un dizionario/lista che gli passo con dentro le informazioni riguardo il tipo di simulazione
           # device = IBMQ.get_provider().get_backend('ibmqx2')
#            properties = device.properties()
#            coupling_map = device.configuration().coupling_map
#            gate_times = [('u1', None, 0), ('u2', None, 60), ('u3', None, 120)]
#            basis_gates = noise_model.basis_gates
#    # Select the QasmSimulator from the Aer provider
#            simulator = Aer.get_backend('qasm_simulator')
            pass
        
            
        def execute(self):
            if self.provider == Aer.get_backend("statevector_simulator"):
                result = execute(self.qcircuit, self.provider).result()
                statevector = result.get_statevector(self.qcircuit)
            elif self.provider == Aer.get_backend("qasm_simulator"):
                result = execute(self.qcircuit, self.IBMInterface.provider, self.IBMInterface.shots, self.IBMInterface.noise_model, self.IBMInterface.noise_model.basis_gates, self.IBMInterface.device.configuration().coupling_map).result()
                counts = result.get_counts(self.qcircuit)
                statevector = af.population_from_counts_dictionary(counts, self.IBMInterface.shots, 2**(len(self.qbits)))
            else:
                result = execute(self.qcircuit, self.IBMInterface.device, shots=self.IBMInterface.shots).result()
                counts = result.get_counts(self.qcircuit)
                statevector = af.population_from_counts_dictionary(counts, self.IBMInterface.shots, 2**(len(self.qbits)))
            return statevector
        
        
class IBMInterface():
    def __init__(self):
        self.provider = None
        self.device = None
        self.shots = None
        self.noise = None
        self.noise_model = None
        
        
    def set_noise_model(self):
        if self.noise:
            self.noise_model = NoiseModel.from_backend(self.device)
            
    def set_device(self, device_string):
        pass
            
#    def set_provider(self, string_provider):
#        self.provider = Aer.get_backend(string_provider)            
        
    

