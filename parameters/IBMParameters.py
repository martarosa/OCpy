#MRqiskit from qiskit.providers.aer.noise import NoiseModel
#MRqiskit from qiskit import IBMQ

class IBMParameters():
    def __init__(self):
        self.quantum_prop_keyword = None
        self.provider = None
        self.device = None
        self.shots = None
        self.noise = None
        self.noise_model = None
        self.coupling_map = None
        
        
    def set_noise_model(self):
        pass
        '''
        #MRqiskit
        if self.noise:
            self.noise_model = NoiseModel.from_backend(self.device.properties())
            self.coupling_map = self.device.configuration().coupling_map
        '''    
    def set_device(self):
        pass
        '''
        #MRqiskit
        provider = IBMQ.load_account()
        self.device = provider.get_backend(self.device)
        '''
