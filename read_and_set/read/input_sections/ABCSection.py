from copy import deepcopy
from abc import ABCMeta
from read_and_set.read import auxiliary_functions as af


#Structure to store and check input files sections.
#Child classes are initialized with
# section name,
# default dictionary,
# allowed values for multiple choice keys
# list of case unsensitive keys


class ABCSection(metaclass=ABCMeta):
    def __init__(self):
        self.section = None
        self.section_default_dictionary = {}
        self.section_dictionary = {}
        self.allowed_val = None
        self.case_unsensitive_keys = None


    def init(self, user_input):
        self.section_dictionary = dict(user_input[self.section])

    def check(self):
        self.check_missing()
        self.lowering_case()
        self.check_allowed_values()

    def fill_empty_with_default(self):
        tmp_dict = deepcopy(self.section_default_dictionary)
        tmp_dict.update(self.section_dictionary)
        self.section_dictionary = tmp_dict

    def check_missing(self):
        for key in self.section_dictionary:
            if key == "missing":
                af.exit_error("ERROR: in namelist " + self.section + " value of " + key + " is missing")


    def check_keys(self):
        if (not all(elem in self.section_dictionary.keys() for elem in self.section_default_dictionary.keys())
                or len(self.section_dictionary.keys()) > len(self.section_default_dictionary.keys())):
            af.exit_error("ERROR in input. Wrong key in " + self.section + "\n")

    def lowering_case(self):
        for key in self.case_unsensitive_keys:
            if key in self.section_dictionary:
                self.section_dictionary[key] = self.section_dictionary[key].lower()


    def check_allowed_values(self):
        if self.allowed_val:
                for key, par in self.allowed_val:
                    if key in self.section_dictionary:
                        if self.section_dictionary[key] not in par:
                            af.exit_error("ERROR in input. Wrong \"" + key + "\" value. \n"
                             + "Default value is: " + self.section_default_dictionary[key] + "\n")


    def check_namelist_key_and_print(self, key, output_string):
        if key in self.section_dictionary:
            print(output_string)


    def check_namelist_key_exist_and_value(self, key, value):
        if key in self.section_dictionary:
            if self.section_dictionary[key] == value:
                return True


    def check_namelist_key_exist_and_list_value(self, key, value):
        if key in self.section_dictionary:
            for i in range(len(value)):
                if self.section_dictionary[key] == value[i]:
                    return True


    def check_namelist_key_exist(self, key):
        if key in self.section_dictionary:
            return True

    def add_key_and_value_to_namelist(self, key, value):
        self.section_dictionary[key] = value


