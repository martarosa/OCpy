from read_and_set.read.input_sections.ABCSection import ABCSection


class SectionRabitzParameters(ABCSection):
    def __init__(self):
        super().__init__()
        self.section = "PARAMETERS"
        self.section_default_dictionary = {'alpha0': '1'}
        self.section_dictionary = {}



