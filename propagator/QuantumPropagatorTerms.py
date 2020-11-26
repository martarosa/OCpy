# -*- coding: utf-8 -*-
"""
Created on Tue Sep 22 11:03:33 2020

@author: castd
"""

import numpy as np

from read_and_set.read import auxiliary_functions as af
from molecule.Molecule import Molecule
#from propagator import math_functions as mf

#MRqiskit from qiskit import QuantumCircuit , QuantumRegister , ClassicalRegister
#MRqiskit from qiskit import Aer , execute , IBMQ
from propagator.ABCPropagatorTerms import ABCPropagatorTerms
#MRqiskit from qiskit.providers.aer.noise import NoiseModel
from parameters.IBMParameters import IBMParameters

import sys


class QuantumPropagatorTerms(ABCPropagatorTerms):
        def __init__(self):
            self.dict_terms = {}
            self.qbits = None
            self.cbits = None
            self.qcircuit = None
            self.IBMParameters = IBMParameters() 

            
        def set_qbits(self, mol):
            self.qbits = QuantumRegister(mol.wf.n_ci, 'q')
            
        def set_cbits(self, mol):
            self.cbits = ClassicalRegister(mol.wf.n_ci, 'c')
            
        def set_qcircuit(self):
            self.qcircuit = QuantumCircuit(self.qbits,self.cbits)
            
            
        def set_qprocessor(self, mol):
            self.set_qbits(mol)
            self.set_cbits(mol)
            self.set_qcircuit()

        def init(self, prop_conf):
            self.IBMParameters.quantum_prop_keyword = prop_conf.quantum_prop_keyword
            self.IBMParameters.provider = prop_conf.provider
            self.IBMParameters.device = prop_conf.device
            self.IBMParameters.shots = prop_conf.shots
            self.IBMParameters.noise = prop_conf.noise
            self.IBMParameters.noise_model = self.IBMParameters.set_noise_model()
            self.dict_terms["quantum_evo"] = self.quantum_evo_circuit
            self.dict_terms["expectation_value"] = self.expectation_value_hermitian_operator
            self.dict_terms["measure_computational_basis"] = self.computational_basis_measurement
            self.dict_terms["prepare_GS_linear_mapping"] = self.prepare_GS_linear_mapping
                    
            
            
##### Method to build the circuit equivalent to simulation of a single_particle hamiltonian, in the field-free diagonal basis, under a TD external perturbation (linear mapping) ##### 

#### N.B. all this section should be updated aware of the new qiskit functionality: WeightedPauliOperator.evolve ######


        def quantum_evo_circuit(self, mol, dt, field_dt_vector):
            for k in range(mol.wf.n_ci):
                self.qcircuit.u1((mol.par.en_ci[k] - np.dot(mol.par.muT,(field_dt_vector)/2)[k][k])*dt, self.qbits[k])
                self.qcircuit.h(self.qbits[k])
            self.bounded_occupation_gatex(mol, dt, field_dt_vector)
            for k in range(mol.wf.n_ci):
                self.qcircuit.rx(np.pi/2, self.qbits[k])
            self.bounded_occupation_gatey(mol, dt, field_dt_vector)
            
        def bounded_occupation_gatex(self, mol, dt, field_dt_vector):
            for i in range(mol.wf.n_ci-1):
                for j in range(1,mol.wf.n_ci):
                    if j > i:
                        self.qcircuit.cx(self.qbits[j], self.qbits[i])
                        self.qcircuit.rz(np.dot(mol.par.muT,((field_dt_vector)/2))[j][i]*dt, self.qbits[i])
                        self.qcircuit.cx(self.qbits[j], self.qbits[i])
                self.qcircuit.h(self.qbits[i])
            self.qcircuit.h(self.qbits[j])
            
        def bounded_occupation_gatey(self, mol, dt, field_dt_vector):
            for i in range(mol.wf.n_ci-1):
                for j in range(1,mol.wf.n_ci):
                    if j > i:
                        self.qcircuit.cx(self.qbits[j], self.qbits[i])
                        self.qcircuit.rz(np.dot(mol.par.muT,((field_dt_vector)/2))[j][i]*dt, self.qbits[i])
                        self.qcircuit.cx(self.qbits[j], self.qbits[i])
                self.qcircuit.rx(-np.pi/2, self.qbits[i])
            self.qcircuit.rx(-np.pi/2, self.qbits[j])

##### Method to build the circuit equivalent to the simulation of a single particle hamiltonian, under a TD external perturbation (binary_mapping) #####

        def binary_quantum_evo(self, mol, dt, perturbation):
            pass	
            
            
###### Measurement_methods #######
                        
        def computational_basis_measurement(self):
            self.qcircuit.measure(self.qbits, self.cbits)
            
        def pure_state_expectation_value(self):
            """ Here, according to the chosen mapping (system --> computer), we compute 
            the projection on an arbitrary pure state as the expectation value of the corresponding 
            density matrix. USE AQUA METHODS """
            pass
            
        def expectation_value_hamiltonian_operator(self):
            """ è verosimile che questo metodo avrà dentro l'inizializzazione di vari circuiti con questo schema:
                CIRCUITO DI EVOLUZIONE + STRINGHE DI PAULI VARIE + MISURA BASE COMPUTAZIONALE
                fino ad esaurimento delle stringhe di pauli.
                Visto che execute è separato deve restituire una lista di circuiti da eseguire
                Per ottenere il valore di aspettazione poi avremo una cosa del tipo:
                for counts in data:
                    expectation_value += counts*matrix_element_giusto """
            ### un'altra possibilità è ricostruire l'hamiltoniano con gli oggetti Pauli di Aqua e utilizzare
            # le funzioni qiskit di expectation_value
            pass
                     
###### Initialization_methods ######
            
        def prepare_GS_linear_mapping(self):
            self.qcircuit.x(self.qbits[0])
            
###### different kinds of execution are available #######
            
        def configure_backend(self, simulator_options):  ## simulator_options è un dizionario/lista che gli passo con dentro le informazioni riguardo il tipo di simulazione
            pass
        
            
        def execute(self):
            print(self.IBMParameters.provider)
            '''
            #MRqiskit
            if self.IBMParameters.provider == "statevector_simulator":
                result = execute(self.qcircuit, Aer.get_backend(self.IBMParameters.provider)).result()
                statevector = result.get_statevector(self.qcircuit)
            elif self.IBMParameters.provider == "qasm_simulator":
                result = execute(self.qcircuit, Aer.get_backend(self.IBMParameters.provider), self.IBMParameters.shots, self.IBMParameters.noise_model, self.IBMParameters.noise_model.basis_gates, self.IBMParameters.device.configuration().coupling_map).result()
                counts = result.get_counts(self.qcircuit)
                statevector = af.population_from_counts_dictionary(counts, self.IBMParameters.shots, 2**(len(self.qbits)))
            else:
                result = execute(self.qcircuit, self.IBMParameters.device, shots=self.IBMParameters.shots).result()
                counts = result.get_counts(self.qcircuit)
                statevector = af.population_from_counts_dictionary(counts, self.IBMParameters.shots, 2**(len(self.qbits)))
            return statevector
            '''
        
          
        
    

