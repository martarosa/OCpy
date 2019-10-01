import numpy as np

from save.SaveOC import SaveOC
from save.SaveTools import SaveTools


class SaveOCGenetic(SaveOC):
    def __init__(self):
        super().__init__()
        self.folder = None
        self.filename = None
        self.restart_calculation = None
        self.save_files = None
        self.restart_file = None
        self.save_tools = SaveTools()


    def init_save_file_list(self, oc_iterator, restart_step):
        #log = SaveFile(".log",
        #               "#Final states populations \n #fields: n_iteration, convergence, J, projection on target state, alpha*integral(field^2)\n",
        #               1,
        #               'log_file',
        #               oc_iterator)

        #self.save_files = [log, final_pop, pop_t, field_t]
        #self.restart_file = SaveRestart("_field_bkp.dat", restart_step, oc_iterator)
        pass



    def init_save(self, save_parameters, log_header_parameters, oc_iterator):
        self.folder = save_parameters.folder
        self.filename = save_parameters.name
        self.restart_calculation = save_parameters.restart_calculation
        self.init_save_file_list(oc_iterator, save_parameters.restart_step)

        for i in np.arange(len(self.save_files)):
            self.save_tools.creation_save_files(self.folder + self.filename, self.save_files[i], self.restart_calculation)
        self.save_tools.print_log_header(self.folder + self.filename, log_header_parameters)


    def save(self, iteration):
        for i in np.arange(len(self.save_files)):
            self.save_tools.save_every_n_iterations(iteration, self.folder + self.filename, self.save_files[i])
        self.save_tools.save_restart(iteration, self.folder + self.filename, self.restart_file)