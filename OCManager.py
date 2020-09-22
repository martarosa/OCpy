import numpy as np

import dictionaries.OCDictionaries
from dictionaries import SaveDictionaries as dict

from SystemObj import Func_tMatrix
from parameters.OCManagerParameters import OCManagerParameters


# OCmanager is only one, and deals with the different OC_iterators (propagation, rapitz, genetic...) Depending
# on the iterator, the SaveOC class is differently initialized. The restart method is inside Save and is also differently initialized
# depending on the algorithm (e.g. for the genetic it reads the field parameters, while for Rabitz it reads the field point by point


class OCManager:
    def __init__(self):
        self.par = OCManagerParameters()

        self.oc_iterator = None
        self.save = None

        self.psi_coeff_t_matrix = Func_tMatrix()
        self.field_psi_matrix = Func_tMatrix()


    def init_oc(self, oc_input, oc_conf, save_input, log_header_input, molecule, starting_field, medium):
        self.par.alpha = oc_input.alpha
        self.par.oc_iterator_name = oc_input.oc_iterator_name
        self.par.convergence_thr = oc_input.convergence_thr
        self.par.n_iterations = oc_input.n_iterations


        self.init_oc_iterator(oc_input,
                              oc_conf,
                              molecule,
                              starting_field,
                              medium,
                              self.set_alpha_t(oc_input.dt, oc_input.nstep))

        self.init_save(save_input, log_header_input)



    def init_oc_iterator(self, oc_input, oc_conf, molecule, starting_field, medium, alpha_t):
        self.oc_iterator = dictionaries.OCDictionaries.OCAlgorithmDict[self.par.oc_iterator_name]()
        self.oc_iterator.init(molecule, starting_field, medium, alpha_t, oc_input, oc_conf)


    def init_save(self, save_parameters, log_header_parameters):
        self.save = dict.SaveDict[self.par.oc_iterator_name]()
        self.save.init_save(save_parameters, log_header_parameters, self.oc_iterator)


    def set_alpha_t(self, dt, nstep):
        alpha_t = np.zeros([nstep])
        for i in range(nstep):
            if self.par.alpha == "const":
                alpha_t[i] = 1  # const
            elif self.par.alpha == "sin":
                alpha_t[i] = 1 / np.square(np.sin(np.pi * (i + 1) / nstep))  # paper gross
            elif self.par.alpha == "quin":
                alpha_t[i] = 1 / np.exp(
                    -np.power((((i) * dt - 125) / 220), 12))  # paper quinolone
        return alpha_t



    def iterate(self):
        current_iteration = 0
        while (current_iteration <= self.par.n_iterations or self.par.convergence_thr < self.oc_iterator.par.convergence_t):
            self.oc_iterator.iterate(current_iteration)
            self.save.save(current_iteration)
            current_iteration += 1
        self.psi_coeff_t_matrix = self.oc_iterator.psi_coeff_t_matrix
        self.field_psi_matrix = self.oc_iterator.field_psi_matrix
        self.par.convergence_t = self.oc_iterator.par.convergence_t



