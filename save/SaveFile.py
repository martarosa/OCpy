class SaveFile():
    def __init__(self, name, header, save_step, dict_flag, oc_iterator):
        self.name = name
        self.header = header
        self.save_step = save_step
        self.dict_flag = dict_flag
        self.out = oc_iterator.dict_out[self.dict_flag]

