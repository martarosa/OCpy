from InitPar import InitPar
from pcm.PCMParameters import PCMParameters
from read.ReadOutputGaussian import ReadOutputGaussian


class InitPCMPar(InitPar):
    def __init__(self):
        super().__init__()
        self.parameters = PCMParameters()

    def init(self, user_input):
        read_output = ReadOutputGaussian()
        self.parameters.cavity = read_output.read_cavity_tesserae(user_input.sys.par['folder'] + user_input.env.par['name_file_cavity'])
        self.parameters.env = user_input.env.par['env']
        if self.parameters.env == 'sol':
             self.parameters.Qnn_reactionfield = read_output.read_Q_matrix(user_input.sys.par['folder'] +
                                                                           user_input.env.par['name_q_tdplas'])
             self.parameters.Qnn_localfield = read_output.read_Q_matrix(user_input.sys.par['folder'] +
                                                                        user_input.env.par['name_q_local_field'])