class SaveFile():
    def __init__(self, name, header, nstep, dict_flag, oc_iterator):
        self.name = name
        self.header = header
        self.nstep = nstep
        self.dict_flag = dict_flag
        self.out = oc_iterator.oc_iterator_parameters.dict_out[self.dict_flag]

