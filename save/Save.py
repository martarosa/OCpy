from save.SaveTools import SaveTools
from abc import ABCMeta, abstractmethod


class Save(metaclass=ABCMeta):

    def __init__(self):
        self.folder = None
        self.filename = None
        self.restart_calculation = None
        self.save_files = None
        self.restart_file = None
        self.save_tools = SaveTools()

    @abstractmethod
    def init_save_file_list(self, oc_iterator, restart_step):
        pass

    @abstractmethod
    def save(self, iteration):
        pass

    @abstractmethod
    def init_save(self, save_input, log_header_input, oc_iterator):
        pass




