# -*- coding: utf-8 -*-

import numpy as np

class ReadOutputPsi4():
    
    def read_psi4_output_file(self, name_file):
        output_file_list_format = np.load(name_file, allow_pickle = True)
        return output_file_list_format
