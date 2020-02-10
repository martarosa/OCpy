from abc import ABCMeta, abstractmethod


class ABCSetInput(metaclass=ABCMeta):
    def __init__(self):
        self.input_parameters = None

    @abstractmethod
    def set(self, user_input):
        pass