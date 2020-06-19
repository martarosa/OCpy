#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  6 15:55:38 2020

@author: davide
"""

import time

import numpy as np
import multiprocessing as mp
import concurrent.futures

from parameters.FieldParameters import FieldParameters
from read_and_set.read import auxiliary_functions as af
from read_and_set.read import NamelistTools as nmtool
from OC.ABCOCIterator import ABCOCIterator
from parameters.OCIteratorParameters import OCIteratorParameters
from SystemObj import DiscreteTimePar
from field.Field import Field, Func_tMatrix
from copy import deepcopy
from scipy import optimize
from propagator import PropagatorsEulero as prop
import sys


import random



class OCScipyOptimizeGradFreeIterator(ABCOCIterator):
    def __init__(self):
        super().__init__()
        self.par = OCIteratorParameters()
        self.field_par = FieldParameters()
        self.discrete_t_par = DiscreteTimePar()
        self.field = Field()
        self.psi_coeff_t_matrix = Func_tMatrix()
        self.dict_out = {}
        self.prop_psi = None
        self.amplitudes_to_save = []
        self.initial_c0 = None
        self.current_iteration = None
        self.result = None
        self.pop_tgt = ['NaN', ]
        
       
        
    def iterate(self, current_iteration):
        if current_iteration == 0:
     #       print(current_iteration)
            self.current_iteration = current_iteration
            self.optimize()  
        else:
            pass
        
   
    
    def initial_guess_control_parameters(self):
        x = np.zeros(self.field_par.fi.size).tolist()
        for k in range(2,5):
            x[k] = 0.01
        return x
    
    def optimize(self):
        result = optimize.minimize(self.calc_J, self.initial_guess_control_parameters(), method = self.par.optimization_method, options = {'disp' : True, 'return_all' : True, 'maxiter' : self.par.n_iterations}, callback = self.callback_for_scipy)                       
        self.result = result                      

        
    def init_output_dictionary(self):
        self.dict_out['log_file'] = self.get_log_file_out
        self.dict_out['field_ampl'] = self.get_field_ampl
        self.dict_out['scipy_result'] = self.result
        
    
    def init(self, molecule, starting_field, pcm, alpha_t, oc_input):
        self.discrete_t_par.nstep = oc_input.nstep
        self.discrete_t_par.dt = oc_input.dt
        self.par.target_state = oc_input.target_state
        self.par.alpha_t = alpha_t
        self.par.convergence_t = 99999
        self.par.J = [99999,]
        self.par.oc_iterator_prop = oc_input.oc_iterator_name
        self.par.propagation_type = oc_input.propagation_type
        self.par.n_iterations = oc_input.n_iterations
        self.par.optimization_method = oc_input.optimization_method

        self.field = deepcopy(starting_field)
        self.field_par = starting_field.par
        self.init_symplex(molecule, starting_field, pcm)
        self.init_output_dictionary()
    
    def check_convergence(self):
        if self.current_iteration - self.par.n_iterations + 1 == 0:
            self.par.convergence_t = self.par.J[0]
        else:
            self.par.convergence_t =  self.par.J[self.current_iteration - self.par.n_iterations + 1] - self.par.J[self.current_iteration - self.par.n_iterations]
            
    
    def calc_J(self, coefficients):
        self.field.par.fi = np.asarray(coefficients).reshape((-1, 3))
        self.field.chose_field('sum', discrete_t_par = self.discrete_t_par)
        if self.par.propagation_type == 'quantum':
           from qiskit.tools.qi.qi import partial_trace
           self.prop_psi.set_propagator(self.prop_psi.propagator_terms.mol, self.prop_psi.propagator_terms.pcm)
           statevector = self.prop_psi.propagate_n_step(self.discrete_t_par, self.field.field)
           p_tgt = np.real(partial_trace(statevector, np.delete(np.arange(len(self.initial_c0)), np.arange(len(self.initial_c0))[np.argmax(self.par.target_state)]))[1,1])
           field = np.real(self.alpha_field_J_integral(self.field.field))
           J = 1 - p_tgt + field
        else:
            self.prop_psi.set_propagator(self.prop_psi.propagator_terms.mol, self.prop_psi.propagator_terms.pcm)
            self.prop_psi.propagator_terms.mol.wf.set_wf(self.initial_c0, 1)
            self.prop_psi.propagate_n_step(self.discrete_t_par, self.field.field)
                       #calcolo popolazione e integrale del campo, per debug
            pop= np.real(af.projector_mean_value(self.prop_psi.propagator_terms.mol.wf.ci,
                                                 self.par.target_state))
            field = np.real(self.alpha_field_J_integral(self.field.field))
           #calcolo J
            J = 1 - np.real(af.projector_mean_value(self.prop_psi.propagator_terms.mol.wf.ci,
                                                self.par.target_state)
                        + self.alpha_field_J_integral(self.field.field))
        return J
        
        
    def chose_propagator(self):
        if self.par.propagation_type == 'quantum':
            from propagator import PropagatorQuantum as qprop
            self.prop_psi = qprop.PropagatorQuantum()
            from qiskit import Aer
            self.prop_psi.provider = Aer.get_backend("statevector_simulator")
        else:
            self.prop_psi = prop.PropagatorEulero2Order()
            
    def init_symplex(self, molecule, starting_field, pcm):
        self.chose_propagator()
        self.prop_psi.set_propagator(molecule, pcm)
        self.initial_c0 = self.prop_psi.propagator_terms.mol.wf.ci
    
    def callback_for_scipy(self, coefficients):
        print(self.current_iteration)
        self.current_iteration += 1
        self.par.J.append(self.calc_J(coefficients))
        self.amplitudes_to_save.append(coefficients)
        print("J at iteration " + str(self.current_iteration) + " is " + str(self.par.J[-1]))
        if self.par.propagation_type == 'classic':
            self.pop_tgt.append(np.real(af.projector_mean_value(self.prop_psi.propagator_terms.mol.wf.ci, self.par.target_state)))
        
            
    def alpha_field_J_integral(self, field):
        ax_square= field.f_xyz.ndim - 1
        ax_integral= field.f_xyz.ndim - 2
        f_square = np.sum(field.f_xyz * field.f_xyz, axis=ax_square) * self.par.alpha_t
        f_integral = np.sum(f_square, axis=ax_integral)
        out_integral = f_integral*self.discrete_t_par.dt
        return out_integral
    

### methods to save the output #####
        
    def get_field_ampl(self):
     #   print(self.current_iteration)
        return np.reshape(self.amplitudes_to_save[self.par.n_iterations - self.current_iteration], ((-1,3)))
    
    
    def field_J_integral(self, amplitudes):
        self.field.par.fi = np.asarray(amplitudes).reshape((-1, 3))
        self.field.chose_field('sum', discrete_t_par = self.discrete_t_par)        
        ax_square = self.field.field.f_xyz.ndim - 1
        ax_integral= self.field.field.f_xyz.ndim - 2
        f_square = np.sum(self.field.field.f_xyz * self.field.field.f_xyz, axis=ax_square)
        f_integral = np.sum(f_square, axis=ax_integral)
        out_field = f_integral*self.discrete_t_par.dt
        return out_field

    
    def get_log_file_out(self): 
        if self.par.propagation_type != 'quantum':
            pop_tgt = self.pop_tgt[self.current_iteration - self.par.n_iterations]
            J = self.par.J[-self.current_iteration]
            integral_field = pop_tgt - J
            self.check_convergence()
            integral_field = self.field_J_integral(self.amplitudes_to_save[self.par.n_iterations - self.current_iteration]) 
            log_array = np.array([self.par.convergence_t, self.par.J[self.par.n_iterations - self.current_iteration], pop_tgt, integral_field])
            self.current_iteration += 1
            print(self.par.convergence_t)
            return log_array
        else:
            self.check_convergence()
            integral_field = self.field_J_integral(self.amplitudes_to_save[self.par.n_iterations - self.current_iteration]) 
            log_array = np.array([self.par.convergence_t, self.par.J[self.current_iteration - self.par.n_iterations + 1], integral_field])
            self.current_iteration += 1
            return log_array

    def get_scipy_result_out(self):
        result = self.result
        return result

    def get_restart(self):
        out = np.concatenate((self.genetic_par.omegas_matrix, np.reshape(self.amplitudes_to_save[-1], ((-1,3)))))
        return out

        
        
