# -*- coding: utf-8 -*-
"""
Created on Fri Nov 27 00:42:51 2020

@author: castd
"""

import numpy as np

import matplotlib.pyplot as plt

from read_and_set.read import auxiliary_functions as af

to_process = np.load("postprocessing_material_bfgs_45par.npy", allow_pickle=True)

to_check = np.load("postprocessing_material_bfgs_45_par_old_implementation.npy", allow_pickle=True)



trajectory = to_process[0]

trajectory_check = to_check[0].f_xyz

field = to_process[1]

field_check = to_check[2]

hamiltonian = to_process[2]

hamiltonian_check = to_check[3]

control_operator_check = to_check[1]

control_operator = to_process[3]

instantaneous_eigenvalues = []

#instantaneous_eigenvalues_check = []

expectation_value_t = []

#expectation_value_t_check = []


for k in range(10000):
    hamiltonian_t =  hamiltonian - control_operator*(1-field[k])
    eigenvalues, eigenvectors = np.linalg.eigh(hamiltonian_t)
  #  eigenvalues_t, eigenvectors_t = np.linalg.eigh(control_operator_check[k])
    instantaneous_eigenvalues.append(eigenvalues)
  #  instantaneous_eigenvalues_check.append(eigenvalues_t)
    expectation_value = af.compute_expectation_value(trajectory[k,:], hamiltonian)
  #  expectation_value_check = af.compute_expectation_value(trajectory_check[k], hamiltonian_check)
    expectation_value_t.append(expectation_value)
  #  expectation_value_t_check.append(expectation_value_check)

