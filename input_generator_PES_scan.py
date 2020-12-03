# -*- coding: utf-8 -*-
"""
Created on Wed Dec  2 09:47:40 2020

@author: castd
"""

import numpy as np
import os
import platform
import shutil


### script to modify and create input.dat files to execute PES scan ###

distance = [distance for distance in np.arange(0.14, 2.14, 0.1)]

path = "/Users/castd/Desktop/OCpy_hybrid_implementation/"  ### modify absolute path at need

input_file_path = path + "input.dat"  #### input file for the calculations to be scan (notice that should be located in the same place of the
                                      #### psi4 output folder)

quantum_chemistry_files = []

for bond_length in distance:
    quantum_chemistry_files.append(path + "psi4output/" + "psi4_output_" + str(bond_length) + "_.npy")

for psi4output in quantum_chemistry_files:
    if platform.system() == 'Windows':
        shutil.copy(input_file_path , input_file_path  + str(distance[quantum_chemistry_files.index(psi4output)]))
    else:
        os.system('cp ' + input_file_path + ' ' + input_file_path + str(distance[quantum_chemistry_files.index(psi4output)]))
    f = open(input_file_path + str(distance[quantum_chemistry_files.index(psi4output)]) , "r+")
    lines = f.readlines()
    for x in lines:
        if x == 'name = \n':
            lines[lines.index(x)] = 'name = H2_scan_BFGS_6_31G_new_implementation_' + str(distance[quantum_chemistry_files.index(psi4output)]) + '\n'
        elif x == 'name_quantum_chemistry = \n':
            lines[lines.index(x)] = 'name_quantum_chemistry =' + ' ' + psi4output + '\n'
    f.close()
    f = open(input_file_path + str(distance[quantum_chemistry_files.index(psi4output)]) , "w")
    f.writelines(lines)
    f.close()
            
        

        



 