from parameters.SaveParameters import SaveParameters
from save.SaveTools import SaveTools
from save.LogHeader import LogHeader
from abc import ABCMeta, abstractmethod


class ABCSave(metaclass=ABCMeta):

    def __init__(self):
        self.par = SaveParameters()

        self.save_files = []                #list of files produced and saved in the case of specific calculation.
                                            # Attributes are assigned upon initialization
        self.restart_file = None
        self.save_tools = SaveTools()

    @abstractmethod
    def init_save_file_list(self, oc_iterator, restart_step):
        pass

    @abstractmethod
    def save(self, iteration):
        pass

    @abstractmethod
    def init_save(self, save_input, log_input, oc_iterator):
        pass


    def print_log_header(self, log_input):
        log_header = LogHeader()
        log_header.init_log_header(log_input)
        f = open(self.par.folder + self.par.filename + ".log", 'a')
        f.write(log_header.header)
        f.close()

