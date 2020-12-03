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
import dictionaries.OCDictionaries as ocdict
import dictionaries.PropagatorDictionaries as pdict
import pandas as pd
import random
import pdb
#MRqiskit from qiskit.quantum_info import partial_trace




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
        self.pop_target = []
        
       
        
    def iterate(self, current_iteration):
        if current_iteration == 0:
            self.current_iteration = current_iteration
            self.optimize()  
        else:
            pass
        
   
    
    def initial_guess_control_parameters(self):
        if self.field.par.field_type == 'gaussian_sum':
            x = np.zeros((int(self.field.par.fi.size/3), 3))
            x[:,1] = np.random.uniform(0, self.discrete_t_par.dt * self.discrete_t_par.nstep, int(self.field.par.fi.size/3))
            x[:,2] = np.random.uniform(-0.01, 0.01, int(self.field.par.fi.size/3))
        elif self.field.par.field_type == "free_harmonics":
            x = np.zeros((int(self.field.par.fi.size/3), 3))
            x[0,2] = -0.9
           # x[:,1] = np.random.uniform(0, 6, int(self.field.par.fi.size/3))
            x[:,2] = np.random.uniform(0, 6, int(self.field.par.fi.size/3))
        x = np.reshape(x, self.field.par.fi.size).tolist()
     #   x = np.loadtxt("/Users/castd/Desktop/Risultati_OC_GS/H2_CISD_space/Genetic_optimization/6-31G/test_new_implementation_genetic_initial_guess_allzero2_gaussian_field_bkp.dat")[10:,:]
        #x = np.reshape(x, self.field.par.fi.size).tolist()
       # x = np.load("/Users/castd/Desktop/Risultati_OC_GS/H2_CISD_space/ScipyOptimizers/BFGS/6-31G/free_sines_perturbation/test_bfgs_single_point_H2_6-31g_free_sines_15par0.74_result.npy", allow_pickle=True)[0].x.tolist()
        return x
    
    ### prendo spunto da restart_genetic (ma meglio), posso dargli
    
    def optimize(self):
        result = optimize.minimize(self.calc_J, self.initial_guess_control_parameters(), method = self.par.oc_iterator_name, options = {'disp' : True, 'return_all' : True, 'maxiter' : self.par.n_iterations, 'gtol': 1e-8}, callback = self.callback_for_scipy)                       
        self.result = result                      

        
    def init_output_dictionary(self):
        self.dict_out['log_file'] = self.get_log_file_out
        self.dict_out['field_ampl'] = self.get_field_ampl
        self.dict_out['field_t'] = self.get_field_t
        self.dict_out['scipy_result'] = self.result
        
    
    def init(self, molecule, starting_field, medium, alpha_t, oc_input, oc_conf, prop_conf):
        self.discrete_t_par.nstep = oc_input.nstep
        self.discrete_t_par.dt = oc_input.dt
        print(self.discrete_t_par.dt)
        self.par.target_state = oc_input.target_state
        self.par.alpha_t = np.array(alpha_t)
        self.par.convergence_t = 99999
        self.par.J = [99999,]
        self.par.control_problem = oc_input.control_problem
        self.obj_fun = ocdict.OCObjectiveFunction[self.par.control_problem]()
        self.obj_fun.set_objective_function()
        self.par.oc_iterator_name = oc_input.oc_iterator_name
        self.par.propagator = oc_input.propagator
        self.par.n_iterations = oc_input.n_iterations

        self.field = deepcopy(starting_field)
        self.field_par = starting_field.par
        self.prop_psi = pdict.PropagatorDict[self.par.propagator]()
        if prop_conf != None:
            prop_conf.quantum_prop_keyword = self.par.propagator
        self.prop_psi.set_propagator(molecule, medium, prop_conf)
        self.init_wf()
        self.init_output_dictionary()
    
    def check_convergence(self):
        if self.current_iteration - self.par.n_iterations + 1 == 0:
          #  print(self.par.n_iterations)
            self.par.convergence_t = self.par.J[0]
        else:
            self.par.convergence_t =  self.par.J[self.current_iteration - self.par.n_iterations] - self.par.J[self.current_iteration - self.par.n_iterations + 1]
            
    
    def calc_J(self, coefficients):
        if self.par.control_problem == "optical_excitation":
            self.field.par.fi = np.asarray(coefficients).reshape((-1, 3))
            self.field.chose_field('sum', discrete_t_par = self.discrete_t_par)
        elif self.par.control_problem == "ground_state":
            if self.field.par.field_type == "simple_sum": ## ?? è così brutto?
                self.field.par.fi = np.asarray(coefficients).reshape((-1, 1))
            else:
                self.field.par.fi = np.asarray(coefficients).reshape((-1, 3))
               # print(self.field.par.fi)
            self.field.chose_field(self.field.par.field_type, discrete_t_par = self.discrete_t_par)
        self.prop_psi.mol.wf.set_wf(self.initial_c0, 1)
        if self.par.propagator == 'quantum_trotter_suzuki':
            pass
            '''
            MRqiskit 
            self.prop_psi.propagator_terms.set_qprocessor(self.prop_psi.mol)
            self.prop_psi.propagate_n_step(self.discrete_t_par, self.field.field)
            if self.prop_psi.propagator_terms.IBMParameters.provider == "statevector_simulator":
               pop_target = np.real(partial_trace(self.prop_psi.final_state_qc, np.delete(np.arange(len(self.initial_c0)), np.arange(len(self.initial_c0))[np.argmax(self.par.target_state)])).data[1,1])
               field = np.real(self.alpha_field_J_integral(self.field.field))
               J = 1 - pop_target + field
            else:
               p_tgt = self.prop_psi.counts_dictionary[np.argmax(self.par.target_state)]
               field = np.real(self.alpha_field_J_integral(self.field.field))
               J = 1 - p_tgt + field
            '''
        else:
            self.prop_psi.propagate_n_step(self.discrete_t_par, self.field.field)
#            trajectory = self.prop_psi.propagate_n_step(self.discrete_t_par, self.field.field)
#            to_save = []
#            to_save.append(trajectory)
#            to_save.append(self.field.field.f_xyz)
#            to_save.append(self.prop_psi.mol.par.hamiltonian)
#            to_save.append(self.prop_psi.mol.par.control_operator)
#            np.save("postprocessing_material_bfgs_45par.npy", to_save)
#            pdb.set_trace()
            # continua qua
            if self.par.control_problem == "optical_excitation":
                self.obj_fun.compute_objective_function(self, self.prop_psi.mol.wf.ci, self.par.target_state)
                J = np.real(self.obj_fun.J)
            elif self.par.control_problem == "ground_state":
                self.obj_fun.compute_objective_function(self, self.prop_psi.mol.wf.ci, self.prop_psi.mol.par.hamiltonian)
                J = np.real(self.obj_fun.J)
        return J
        
            
    def init_wf(self):
        if self.par.control_problem == 'optical_excitation':
            self.initial_c0 = self.prop_psi.mol.wf.ci
        elif self.par.control_problem == 'ground_state':
            self.initial_c0 = np.zeros(len(self.prop_psi.mol.par.basis_determinants))
            self.initial_c0[-1] = 1
    
    def callback_for_scipy(self, coefficients):
        print(self.current_iteration)
        self.current_iteration += 1
        self.par.J.append(self.calc_J(coefficients))
        self.amplitudes_to_save.append(coefficients)
        print("J at iteration " + str(self.current_iteration) + " is " + str(self.par.J[-1]))
        if self.par.propagator != 'quantum_trotter_suzuki':
            pass
    #        self.pop_target.append(np.real(af.projector_mean_value(self.prop_psi.mol.wf.ci, self.par.target_state)))
        elif self.prop_psi.propagator_terms.IBMParameters.provider == "statevector_simulator":
            pass
            #MRqiskit self.pop_tgt.append(np.real(partial_trace(self.prop_psi.final_state_qc, np.delete(np.arange(len(self.initial_c0)), np.arange(len(self.initial_c0))[np.argmax(self.par.target_state)])).data[1,1]))
        else:
            pass
     #       self.pop_tgt.append(self.prop_psi.counts_dictionary[np.argmax(self.par.target_state)])
        if self.current_iteration == self.par.n_iterations - 1:
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
        print(len(self.amplitudes_to_save))
        print(self.current_iteration - self.par.n_iterations + 1)
        if self.current_iteration - self.par.n_iterations + 1 < len(self.amplitudes_to_save):
            return np.reshape(self.amplitudes_to_save[self.current_iteration - self.par.n_iterations + 1], ((-1,3)))
        else:
            return np.reshape(self.amplitudes_to_save[-1], ((-1,3)))
    
    
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
        if self.par.control_problem == 'optical_excitation':
            if self.current_iteration - self.par.n_iterations + 1 < len(self.pop_target):
                pop_target = self.pop_target[self.current_iteration - self.par.n_iterations + 1]
                integral_field = self.field_J_integral(self.amplitudes_to_save[self.current_iteration - self.par.n_iterations + 1]) 
            else:
                pop_target = self.pop_target[-1]
                integral_field = self.field_J_integral(self.amplitudes_to_save[-1])
            self.check_convergence()
            log_array = np.array([self.par.convergence_t, self.par.J[self.current_iteration - self.par.n_iterations + 1], pop_target, integral_field])
        elif self.par.control_problem == 'ground_state':
            self.check_convergence()
            log_array = np.array([self.par.convergence_t, self.par.J[self.current_iteration - self.par.n_iterations + 1]])
        self.current_iteration += 1
        return log_array


    def get_scipy_result_out(self):
        result = self.result
        return result
    
    def get_restart(self):
        separate1 = np.array([["###","omegas","###"],["###","###","###"]])
        separate2 = np.array([["###", "fi", "###"],["###","###","###"]])
        out = pd.DataFrame(separate1)
        out_tmp = pd.DataFrame(self.field.par.omega)
        out = out.append(out_tmp)
        out_tmp = pd.DataFrame(separate2)
        out = out.append(out_tmp)
        out_tmp = pd.DataFrame(np.reshape(self.amplitudes_to_save[-1], ((-1,3))))
        out = out.append(out_tmp)
        return out

    
    def get_field_t(self):
        print(self.field_psi_matrix.f_xyz)
        field_t_matrix = np.insert(self.field_psi_matrix.f_xyz, 0, self.field_psi_matrix.time_axis, axis = 1)
        return field_t_matrix

        
        
