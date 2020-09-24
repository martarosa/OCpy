import dictionaries.OCDictionaries
from alpha.Alpha import Alpha
from dictionaries import SaveDictionaries as sdict

from parameters.OCManagerParameters import OCManagerParameters


# OCmanager has only one instance, store all the data and performs all the calculations.
# In particular it initialize the OCIterator() and the SaveOC()
# Depending on the iterator, the SaveOC class is differently initialized.
# Methods for restart are inside SaveOC


class OCManager:
    def __init__(self):
        self.par = OCManagerParameters()

        self.oc_iterator = None
        self.save = None


    def init_oc(self, oc_input, oc_conf, alpha_input, save_input, log_header_input, molecule, starting_field, medium):
        self.par.alpha = oc_input.alpha
        self.par.oc_iterator_name = oc_input.oc_iterator_name
        self.par.convergence_thr = oc_input.convergence_thr
        self.par.n_iterations = oc_input.n_iterations

        alpha = Alpha()
        alpha.init_alpha(alpha_input)
        self.init_oc_iterator(oc_input,
                              oc_conf,
                              molecule,
                              starting_field,
                              medium,
                              alpha.alpha_t)

        self.init_save(save_input, log_header_input)



    def init_oc_iterator(self, oc_input, oc_conf, molecule, starting_field, medium, alpha_t):
        self.oc_iterator = dictionaries.OCDictionaries.OCAlgorithmDict[self.par.oc_iterator_name]()
        self.oc_iterator.init(molecule, starting_field, medium, alpha_t, oc_input, oc_conf)


    def init_save(self, save_parameters, log_header_parameters):
        self.save = sdict.SaveDict[self.par.oc_iterator_name]()
        self.save.init_save(save_parameters, log_header_parameters, self.oc_iterator)


    def iterate(self):
        current_iteration = 0
        while (current_iteration <= self.par.n_iterations or self.par.convergence_thr < self.oc_iterator.par.convergence_t):
            self.oc_iterator.iterate(current_iteration)
            self.save.save(current_iteration)
            current_iteration += 1


