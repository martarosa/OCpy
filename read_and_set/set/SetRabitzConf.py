from read_and_set.set.ABCSetInput import ABCSetInput
from read_and_set.input.OCRabitzInput import OCRabitzInput


class SetRabitzOCInput(ABCSetInput):
    def __init__(self):
        super().__init__()
        self.input_parameters = OCRabitzInput()

    def set(self, user_input):
        self.input_parameters.alpha0 = int(user_input.parameters.section_dictionary["alpha0"])
