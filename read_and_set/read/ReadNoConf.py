from read_and_set.read.ABCReadInputFile import ABCReadInputFile

class ReadNoConf(ABCReadInputFile):
    def __init__(self):
        super().__init__()



    def read_file(self, folder, namefile):
        pass


    def set_nml_sections(self, user_input):
        pass