import numpy as np

from save.SaveFile import SaveFile
from save.ABCSave import ABCSave, SaveParameters
from save.SaveRestart import SaveRestart
from save.SaveTools import SaveTools


class SaveEulero(ABCSave):
    def __init__(self):
        super().__init__()
        self.par = SaveParameters()

        self.save_files = None
        self.restart_file = None
        self.save_tools = SaveTools()

    def init_save_file_list(self, oc_iterator, restart_step):
        pop_t = SaveFile("_pop_t.dat",
                         "#Final states populations \n #fields: n_iteration, nstep, time, states population  \n",
                         restart_step,
                         'pop_t',
                         oc_iterator)
        self.save_files = [pop_t]
        self.restart_file = SaveRestart("_field_bkp.dat", restart_step, oc_iterator)


    def init_save(self, save_input, log_input, oc_iterator):
        self.par.folder = save_input.folder
        self.par.filename = save_input.name
        self.par.restart_calculation = save_input.restart_calculation
        self.par.append = save_input.append
        self.init_save_file_list(oc_iterator, save_input.restart_step)
        for i in np.arange(len(self.save_files)):
            self.save_tools.creation_save_files(self.par.folder + self.par.filename, self.save_files[i], self.par.append)



    def save(self, iteration):
        for i in np.arange(len(self.save_files)):
            self.save_tools.save_every_n_iterations(iteration, self.par.folder + self.par.filename, self.save_files[i])
        self.save_tools.save_restart(iteration, self.par.folder + self.par.filename, self.restart_file)