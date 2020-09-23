from read_and_set.set.ABCSetInput import ABCSetInput


class SetNoConf(ABCSetInput):
    def __init__(self):
        super().__init__()
        self.input_parameters = None


    def set(self, user_input):
        pass