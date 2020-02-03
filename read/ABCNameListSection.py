import sys
from abc import ABCMeta


class ABCNamelistSection(metaclass=ABCMeta):
    def __init__(self):
        self.section = None
        self.section_default_dictionary = {}
        self.section_dictionary = {}
        self.allowed_val = None
        self.case_unsensitive_keys = None


    def init(self, user_input, default, section):
        self.section_dictionary = dict(user_input[section])
        self.section_default_dictionary = default[section]


    def check(self):
        self.check_keys()
        self.lowering_case()
        self.check_allowed_values()

    def check_keys(self):
            if not all(elem in self.section_dictionary.keys() for elem in self.section_default_dictionary.keys()):
                sys.exit("Error in input. Wrong key in " + self.section + "\n")

    def lowering_case(self):
        for key in self.case_unsensitive_keys:
            self.section_dictionary[key] = self.section_dictionary[key].lower()

    def check_allowed_values(self):
        if self.allowed_val:
            for key, par in self.allowed_val:
                if self.section_dictionary[key] not in par:
                    sys.exit("Error in input. Wrong \"" + key + "\" value. \n"
                             + "Default value is: " + self.section_default_dictionary[key] + "\n")