import numpy as np

class MoleculeInput:
    def __init__(self):
        self.wf_ci = None
        self.muT = None
        self.en_ci = None
        self.Vijn = np.array([None])