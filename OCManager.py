import numpy as np

from OCIterator import OCIterator
from OCEuleroIterator import  Eulero1PropagationIterator, Eulero2PropagationIterator
from OCRabitzIterator import OCRabitzIterator
#from OCGeneticIterator import OCGeneticIterator
from save.SaveOC import SaveOC
from save.SaveOCRabitz import SaveOCRabitz
from save.SaveEulero import SaveEulero
from save.SaveOCGenetic import SaveOCGenetic
from OCGeneticIterator import OCGeneticIterator


class OCManager:
    def __init__(self):

        self.oc_iterator_name = None
        self.convergence_thr = None
        self.n_iterations = None

        self.alpha_0 = None
        self.alpha = None

        self.oc_iterator = OCIterator()
        self.save = SaveOC()

        self.psi_coeff_t = None
        self.field_psi_matrix = None
        self.current_iteration = 0
        self.convergence_t = 99999

        #self.restart = None


    def init_oc(self, oc_parameters, iterator_parameters, save_parameters, log_header_parameters, molecule, starting_field, pcm):
        self.alpha0 = oc_parameters.alpha0
        self.alpha = oc_parameters.alpha
        self.oc_iterator_name = oc_parameters.oc_iterator_name
        if self.oc_iterator_name != "eulero_1order" and self.oc_iterator_name != "eulero_2order":
            self.convergence_thr = oc_parameters.convergence_thr
            self.n_iterations = oc_parameters.n_iterations
        else:
            self.n_iterations = 0
            self.convergence_thr = 99999

        self.init_oc_iterator(oc_parameters,
                              iterator_parameters,
                              molecule,
                              starting_field,
                              pcm,
                              self.set_alpha_t(oc_parameters.nstep, oc_parameters.dt))

        self.init_save(save_parameters, log_header_parameters)
        self.field_psi_matrix = np.copy(starting_field.field)



    def init_oc_iterator(self, oc_parameters, iterator_parameters, molecule, starting_field, pcm, alpha_t):
        if self.oc_iterator_name == "rabitzi" or self.oc_iterator_name == "rabitzii":
            self.oc_iterator = OCRabitzIterator()
        elif self.oc_iterator_name == "genetic":
            self.oc_iterator = OCGeneticIterator()
        elif self.oc_iterator_name == "eulero_1order":
            self.oc_iterator = Eulero1PropagationIterator()
        elif self.oc_iterator_name == "eulero_2order":
            self.oc_iterator = Eulero2PropagationIterator()
        self.oc_iterator.init(oc_parameters, iterator_parameters, molecule, starting_field, pcm, alpha_t)



    def init_save(self, save_parameters, log_header_parameters):
        if self.oc_iterator_name == "rabitzi" or self.oc_iterator_name == "rabitzii":
            self.save = SaveOCRabitz()
            self.save.init_save(save_parameters, log_header_parameters, self.oc_iterator)
        elif self.oc_iterator_name == "eulero_1order" or self.oc_iterator_name == "eulero_2order":
            self.save = SaveEulero()
            self.save.init_save(save_parameters, log_header_parameters, self.oc_iterator)
        elif self.oc_iterator_name == "genetic":
            self.save = SaveOCGenetic()
            self.save.init_save(save_parameters, log_header_parameters, self.oc_iterator)


    def set_alpha_t(self, nstep, dt):
        alpha_t = np.zeros([nstep])
        for i in range(nstep):
            if self.alpha == "const":
                alpha_t[i] = self.alpha0  # const
            elif self.alpha == "sin":
                alpha_t[i] = self.alpha0 / np.square(np.sin(np.pi * (i + 1) / nstep))  # paper gross
            elif self.alpha == "quin":
                alpha_t[i] = self.alpha0 / np.exp(
                    -np.power((((i) * dt - 125) / 220), 12))  # paper quinolone
        return alpha_t




    def iterate(self):
        while (self.current_iteration <= self.n_iterations or self.convergence_thr < self.convergence_t):
            self.oc_iterator.iterate(self.current_iteration)
            self.save.save(self.current_iteration)
            self.convergence_t = self.oc_iterator.class_attributes.convergence_t
            self.current_iteration += 1
        self.psi_coeff_t = self.oc_iterator.class_attributes.psi_coeff_t
        self.field_psi_matrix = self.oc_iterator.class_attributes.field_psi_matrix
        self.convergence_t = self.oc_iterator.class_attributes.convergence_t
























