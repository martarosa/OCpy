import numpy as np


from OCEuleroIterator import  Eulero1PropagationIterator, Eulero2PropagationIterator
from OCRabitzIterator import OCRabitzIterator

from save.SaveOCRabitz import SaveOCRabitz
from save.SaveEulero import SaveEulero
from save.SaveOCGenetic import SaveOCGenetic
from OCGeneticIterator import OCGeneticIterator

# OCmanager is only one, and deals with the different OC_iterators (propagation, rapitz, genetic...) Depending
# on the iterator, the SaveOC class is differently initialized. The restart method is inside Save and is also differently initialized
# depending on the algorithm (e.g. for the genetic it reads the field parameters, while for Rabitz it reads the field point by point

class OCManagerParameters:
    def __init__(self):
        self.oc_iterator_name = None
        self.convergence_thr = None
        self.n_iterations = None
        self.alpha0 = None
        self.alpha = None



class OCManager:
    def __init__(self):


 #       self.oc_iterator_name = None
 #       self.convergence_thr = None
 #       self.n_iterations = None
 #       self.alpha_0 = None
 #       self.alpha = None

        self.par = OCManagerParameters()

        self.oc_iterator = None #OCIterator()
        self.save = None #Save()


        self.psi_coeff_t = None
        self.field_psi_matrix = None

        self.current_iteration = 0


    def init_oc(self, oc_input, save_input, log_header_input, molecule, starting_field, pcm):

        self.par.alpha0 = oc_input.alpha0
        self.par.alpha = oc_input.alpha
        self.par.oc_iterator_name = oc_input.oc_iterator_name
        if self.par.oc_iterator_name != "eulero_1order_prop" and self.par.oc_iterator_name != "eulero_2order_prop":
            self.par.convergence_thr = oc_input.convergence_thr
            self.par.n_iterations = oc_input.n_iterations
        else:
            self.par.n_iterations = 0
            self.par.convergence_thr = 99999


        self.init_oc_iterator(oc_input,
                              molecule,
                              starting_field,
                              pcm,
                              self.set_alpha_t(oc_input.dt, oc_input.nstep))

        self.init_save(save_input, log_header_input)
        #self.field_psi_matrix = np.copy(starting_field.get_full_field())


    def init_oc_iterator(self, oc_input, molecule, starting_field, pcm, alpha_t):
        if self.par.oc_iterator_name == "rabitzi" or self.par.oc_iterator_name == "rabitzii":
            self.oc_iterator = OCRabitzIterator()
        elif self.par.oc_iterator_name == "genetic":
            self.par.oc_iterator = OCGeneticIterator()
        elif self.par.oc_iterator_name == "eulero_1order_prop":
            self.par.oc_iterator = Eulero1PropagationIterator()
        elif self.par.oc_iterator_name == "eulero_2order_prop":
            self.par.oc_iterator = Eulero2PropagationIterator()
        self.par.oc_iterator.init(oc_input, molecule, starting_field, pcm, alpha_t)



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
        self.psi_coeff_t = self.oc_iterator.oc_iterator_parameters.psi_coeff_t
        self.field_psi_matrix = self.oc_iterator.oc_iterator_parameters.field_psi_matrix
        self.convergence_t = self.oc_iterator.oc_iterator_parameters.convergence_t













