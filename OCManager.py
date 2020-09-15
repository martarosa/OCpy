import numpy as np


from OC.OCEuleroIterator import  Eulero1PropagationIterator, Eulero2PropagationIterator
from OC.OCRabitzIterator import OCRabitzIterator

from SystemObj import Func_tMatrix
from parameters.OCManagerParameters import OCManagerParameters

from save.SaveOCRabitz import SaveOCRabitz
from save.SaveEulero import SaveEulero
from save.SaveOCGenetic import SaveOCGenetic
from OC.OCGeneticIterator import OCGeneticIterator

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

        self.current_iteration = 0


    def init_oc(self, oc_input, iterator_config_input, save_input, log_header_input, molecule, starting_field, pcm):
        self.par.alpha0 = oc_input.alpha0
        self.par.alpha = oc_input.alpha
        self.par.oc_iterator_name = oc_input.oc_iterator_name
        if self.par.oc_iterator_name != "none":
            self.par.convergence_thr = oc_input.convergence_thr
            self.par.n_iterations = oc_input.n_iterations
        else:
            self.par.n_iterations = 0
            self.par.convergence_thr = 99999


        self.init_oc_iterator(oc_input,
                              iterator_config_input,
                              molecule,
                              starting_field,
                              pcm,
                              self.set_alpha_t(oc_input.dt, oc_input.nstep))

        self.init_save(save_input, log_header_input)






    def init_oc_iterator(self, oc_input, iterator_config_input, molecule, starting_field, pcm, alpha_t):
        if self.par.oc_iterator_name == "rabitzi" or self.par.oc_iterator_name == "rabitzii":
            self.oc_iterator = OCRabitzIterator()
            self.oc_iterator.init(molecule, starting_field, pcm, alpha_t, oc_input)
        elif self.par.oc_iterator_name == "genetic":
            self.oc_iterator = OCGeneticIterator()
            self.oc_iterator.init(molecule, starting_field, pcm, alpha_t, oc_input, iterator_config_input)
        elif self.par.oc_iterator_name == "eulero_1order_prop":
            self.oc_iterator = Eulero1PropagationIterator()
            self.oc_iterator.init(molecule, starting_field, pcm, alpha_t, oc_input)
        elif self.par.oc_iterator_name == "eulero_2order_prop":
            self.oc_iterator = Eulero2PropagationIterator()
            self.oc_iterator.init(molecule, starting_field, pcm, alpha_t, oc_input)



    def init_save(self, save_parameters, log_header_parameters):
        if self.par.oc_iterator_name == "rabitzi" or self.par.oc_iterator_name == "rabitzii":
            self.save = SaveOCRabitz()
            self.save.init_save(save_parameters, log_header_parameters, self.oc_iterator)
        elif self.par.oc_iterator_name == "eulero_1order_prop" or self.par.oc_iterator_name == "eulero_2order_prop":
            self.save = SaveEulero()
            self.save.init_save(save_parameters, log_header_parameters, self.oc_iterator)
        elif self.par.oc_iterator_name == "genetic":
            self.save = SaveOCGenetic()
            self.save.init_save(save_parameters, log_header_parameters, self.oc_iterator)


    def set_alpha_t(self, dt, nstep):
        alpha_t = np.zeros([nstep])
        for i in range(nstep):
            if self.par.alpha == "const":
                alpha_t[i] = self.par.alpha0  # const
            elif self.par.alpha == "sin":
                alpha_t[i] = self.par.alpha0 / np.square(np.sin(np.pi * (i + 1) / nstep))  # paper gross
            elif self.par.alpha == "quin":
                alpha_t[i] = self.par.alpha0 / np.exp(
                    -np.power((((i) * dt - 125) / 220), 12))  # paper quinolone
        return alpha_t


    def iterate(self):
        while (self.current_iteration <= self.par.n_iterations or self.par.convergence_thr < self.oc_iterator.par.convergence_t):
            self.oc_iterator.iterate(self.current_iteration)
            self.save.save(self.current_iteration)
            self.current_iteration += 1
        self.psi_coeff_t_matrix = self.oc_iterator.psi_coeff_t_matrix
        self.field_psi_matrix = self.oc_iterator.field_psi_matrix
        self.par.convergence_t = self.oc_iterator.par.convergence_t













