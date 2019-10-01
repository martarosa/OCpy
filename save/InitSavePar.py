from InitPar import InitPar
from save.SaveParameters import SaveParameters


class InitSavePar(InitPar):
    def __init__(self):
        super().__init__()
        self.parameters = SaveParameters()

    def init(self, user_input):
        self.parameters.folder = user_input.sys.par['folder']
        self.parameters.name = user_input.sys.par['name']
        self.parameters.restart_calculation = user_input.oc.par['restart']
        self.parameters.restart_step = int(user_input.save.par['restart_step'])