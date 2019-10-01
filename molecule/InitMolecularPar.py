from InitPar import InitPar
from molecule.MoleculeParameters import MoleculeParameters
from read.ReadOutputGaussian import ReadOutputGaussian


class InitMolecularPar(InitPar):
    def __init__(self):
        super().__init__()
        self.parameters = MoleculeParameters()


    def init(self, user_input):
        read_output = ReadOutputGaussian()

        self.parameters.wf_ci = read_output.read_ci0(user_input.sys.par['folder'] +
                                                     user_input.wf.par['name_ci'])

        self.parameters.en_ci = read_output.read_en_ci0(user_input.sys.par['folder'] +
                                                        user_input.wf.par['name_ei'])
        self.parameters.muT = read_output.read_muT(user_input.sys.par['folder'] +
                                                   user_input.wf.par['name_mut'], self.parameters.wf_ci.size)
        if user_input.env.par['env'] != 'vac':
            self.parameters.Vijn = read_output.read_V(user_input.sys.par['folder'] +
                                                      user_input.env.par['name_vij'],
                                                      self.parameters.wf_ci.size)