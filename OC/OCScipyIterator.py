#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  6 15:55:38 2020

@author: davide
"""

import time

import numpy as np


from parameters.FieldParameters import FieldParameters
from read_and_set.read import auxiliary_functions as af
from OC.ABCOCIterator import ABCOCIterator
from parameters.OCIteratorParameters import OCIteratorParameters
from SystemObj import DiscreteTimePar
from field.Field import Field, Func_tMatrix
from copy import deepcopy
from scipy import optimize
import dictionaries.PropagatorDictionaries as pdict
from qiskit.quantum_info import partial_trace



import random



class OCScipyOptimizeIterator(ABCOCIterator):
    def __init__(self):
        super().__init__()
        self.par = OCIteratorParameters()
        self.field_par = FieldParameters()
        self.discrete_t_par = DiscreteTimePar()
        self.field = Field()
        self.field_psi_matrix = Func_tMatrix()
        self.dict_out = {}
        self.prop_psi = None
        self.amplitudes_to_save = []
        self.initial_c0 = None
        self.current_iteration = None
        self.result = None
        self.pop_tgt = ['NaN', ]
        
       
        
    def iterate(self, current_iteration):
        if current_iteration == 0:
            self.current_iteration = current_iteration
            self.optimize()  
        else:
            pass
        
   
    
    def initial_guess_control_parameters(self):
#        x = np.zeros(self.field_par.fi.size).tolist()
#        for k in range(15,18):
#            x[k] = 0.01
#        x = [[-3.68766585e-04,  5.33849709e-04,  1.13451009e-06],
# [-3.74400120e-04, -2.16517508e-04, -3.58753124e-04],
# [-1.28134237e-04,  2.05118932e-03, -2.90582851e-05],
# [ 2.96332311e-04, -2.60628692e-04, -7.75910515e-05],
# [-2.90965929e-04, -1.70460776e-03,  9.84011298e-05],
# [-1.57238666e-03,  9.94663608e-03, -4.15938655e-05],
# [-3.65503510e-04, -7.38272642e-05,  3.36194159e-04],
# [-1.15521010e-03, -8.05589509e-04,  1.64288026e-04]]
#        x = [-0.0017128658964207931, 0.001672121823720208, -0.0008831932556318666,
#0.0013729537110224732, 0.0034582972249756036, 0.004595770411588711,
#-0.0025069416233003184, -0.0010888502697676092, 0.003656605855388988,
#0.0011180412233331965, 0.0029888073391043583, -0.003592113597263398,
#0.0040112177903720205, 0.00016570575047176675, -0.0025546752290189507,
#-0.003420106944052337, 0.0027365733243993775, 0.0017089035454578762,
#-0.0038862148527844895, -0.0026, 0.0005755712754066345,
#0.0049, 0.00127372281785325, -0.0008158021974123711]
#        x = [0.0017923175578415933, 0.0014335930579977962, -0.002655151736280887,
#0.002128291843801767, -0.00032958081447407493, 0.0014600087773807564,
#-0.0021757068123361015, -0.005, 0.0006704949806760397,
#-0.005, 0.0049, 0.000371260154670704,
#0.002631845155916659, -0.003465869707698245, -0.0021577813203554545,
#0.0023142461352720165, -0.0006030310331914994, -0.0015803229995185394,
#0.0013317549883324828, 0.003803536060976363, 0.001018553462761253,
#-0.0008046712325593857, -0.0006385955762696089, 0.0032563892888197946]
        x = [-0.004826738588437056, 0.005, -0.0022755529889474657,
-0.0017997225621921528, -0.003674379492403218, -0.0037,
-0.0027, 0.0027950558838909963, 0.003409858220055382,
-0.0037176117298505723, 0.001, -0.0043,
-0.0031284122249297767, 0.005, 0.0009019813759669343,
0.0017008287045830017, 0.001484797880927386, 0.00438139828967898,
0.0029638696380972196, -0.001157443795742901, -0.0022300823916257245,
0.005, -0.003984008037552428, 0.0002927285381301723]
        return x
    
    ### prendo spunto da restart_genetic (ma meglio), posso dargli
    
    def optimize(self):
        result = optimize.minimize(self.calc_J, self.initial_guess_control_parameters(), method = self.par.oc_iterator_name, options = {'disp' : True, 'return_all' : True, 'maxiter' : self.par.n_iterations}, callback = self.callback_for_scipy)                       
        self.result = result                      

        
    def init_output_dictionary(self):
        self.dict_out['log_file'] = self.get_log_file_out
        self.dict_out['field_ampl'] = self.get_field_ampl
        self.dict_out['field_t'] = self.get_field_t
        self.dict_out['scipy_result'] = self.result
        
    
    def init(self, molecule, starting_field, medium, alpha_t, oc_input, oc_conf, prop_conf):
        self.discrete_t_par.nstep = oc_input.nstep
        self.discrete_t_par.dt = oc_input.dt
        self.par.target_state = oc_input.target_state
        self.par.alpha_t = np.array(alpha_t)
        self.par.convergence_t = 99999
        self.par.J = [99999,]
        self.par.oc_iterator_name = oc_input.oc_iterator_name
        self.par.propagator = oc_input.propagator
        self.par.n_iterations = oc_input.n_iterations

        self.field = deepcopy(starting_field)
        self.field_par = starting_field.par
        self.prop_psi = pdict.PropagatorDict[self.par.propagator]()
        prop_conf.quantum_prop_keyword = self.par.propagator
        self.prop_psi.set_propagator(molecule, medium, prop_conf)
        self.init_wf()
        self.init_output_dictionary()
    
    def check_convergence(self):
        if self.current_iteration - self.par.n_iterations + 1 == 0:
            self.par.convergence_t = self.par.J[0]
        else:
            self.par.convergence_t =  self.par.J[self.current_iteration - self.par.n_iterations] - self.par.J[self.current_iteration - self.par.n_iterations + 1]
            
    
    def calc_J(self, coefficients):
        self.field.par.fi = np.asarray(coefficients).reshape((-1, 3))
        self.field.chose_field('sum', discrete_t_par = self.discrete_t_par)
        self.prop_psi.mol.wf.set_wf(self.initial_c0, 1)
        if self.par.propagator == 'quantum_trotter_suzuki':
            self.prop_psi.propagator_terms.set_qprocessor(self.prop_psi.mol)
            self.prop_psi.propagate_n_step(self.discrete_t_par, self.field.field)
            if self.prop_psi.propagator_terms.IBMParameters.provider == "statevector_simulator":
               p_tgt = np.real(partial_trace(self.prop_psi.final_state_qc, np.delete(np.arange(len(self.initial_c0)), np.arange(len(self.initial_c0))[np.argmax(self.par.target_state)])).data[1,1])
               field = np.real(self.alpha_field_J_integral(self.field.field))
               J = 1 - p_tgt + field
            else:
               p_tgt = self.prop_psi.counts_dictionary[np.argmax(self.par.target_state)]
               field = np.real(self.alpha_field_J_integral(self.field.field))
               J = 1 - p_tgt + field
        else:
            self.prop_psi.propagate_n_step(self.discrete_t_par, self.field.field)
            p_tgt = np.real(af.projector_mean_value(self.prop_psi.mol.wf.ci, self.par.target_state))
            field = np.real(self.alpha_field_J_integral(self.field.field))
            J = 1 - p_tgt + field
        return J
        
            
    def init_wf(self):
        self.initial_c0 = self.prop_psi.mol.wf.ci
    
    def callback_for_scipy(self, coefficients):
        print(self.current_iteration)
        self.current_iteration += 1
        self.par.J.append(self.calc_J(coefficients))
        self.amplitudes_to_save.append(coefficients)
        print("J at iteration " + str(self.current_iteration) + " is " + str(self.par.J[-1]))
        if self.par.propagator != 'quantum_trotter_suzuki':
            self.pop_tgt.append(np.real(af.projector_mean_value(self.prop_psi.mol.wf.ci, self.par.target_state)))
        elif self.prop_psi.propagator_terms.IBMParameters.provider == "statevector_simulator":
            self.pop_tgt.append(np.real(partial_trace(self.prop_psi.final_state_qc, np.delete(np.arange(len(self.initial_c0)), np.arange(len(self.initial_c0))[np.argmax(self.par.target_state)])).data[1,1]))
        else:
            self.pop_tgt.append(self.prop_psi.counts_dictionary[np.argmax(self.par.target_state)])
        if self.current_iteration == self.par.n_iterations:
            self.field_psi_matrix = self.field.field
            
            
        
            
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
        if self.par.propagator != 'quantum_trotter_suzuki':
            pop_tgt = self.pop_tgt[self.current_iteration - self.par.n_iterations]
            J = self.par.J[-self.current_iteration]
        #    integral_field = pop_tgt - str(J)
            self.check_convergence()
            integral_field = self.field_J_integral(self.amplitudes_to_save[self.par.n_iterations - self.current_iteration]) 
            log_array = np.array([self.par.convergence_t, self.par.J[self.par.n_iterations - self.current_iteration], pop_tgt, integral_field])
            self.current_iteration += 1
            print(self.par.convergence_t)
            return log_array
        else:
            self.check_convergence()
            integral_field = self.field_J_integral(self.amplitudes_to_save[self.par.n_iterations - self.current_iteration]) 
            log_array = np.array([self.par.convergence_t, self.par.J[self.current_iteration - self.par.n_iterations + 1], self.pop_tgt[self.current_iteration - self.par.n_iterations + 1] ,integral_field])
            self.current_iteration += 1
            return log_array

    def get_scipy_result_out(self):
        result = self.result
        return result

    def get_restart(self):
        out = np.concatenate((self.genetic_par.omegas_matrix, np.reshape(self.amplitudes_to_save[-1], ((-1,3)))))
        return out
    
    def get_field_t(self):
        field_t_matrix = np.insert(self.field_psi_matrix.f_xyz, 0, self.field_psi_matrix.time_axis, axis = 1)
        return field_t_matrix

        
        
