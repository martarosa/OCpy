from qiskit.providers.aer.noise import NoiseModel


class IBMParameters():
    def __init__(self):
        self.quantum_prop_keyword = None
        self.provider = None
        self.device = None
        self.shots = None
        self.noise = None
        self.noise_model = None
        
        
    def set_noise_model(self):
        if self.noise:
            self.noise_model = NoiseModel.from_backend(self.device)
            
    def set_device(self, device_string):
        pass
    
