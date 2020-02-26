import numpy as np


class WaveFunction:
    def __init__(self):
        self.ci = None  # wavefunction coefficients at a given time
        self.ci_prev = None  # past wavefunction coefficients nx2
        self.n_ci = None

    def set_wf(self, ci, norm):
        self.ci = ci
        self.n_ci = ci.size
        if norm:
            self.ci = ci / np.sqrt(np.dot(np.conj(ci), ci))

        self.ci_prev = np.zeros([2, self.n_ci], dtype=complex)
        self.ci_prev[0] = np.copy(self.ci)
        self.ci_prev[1] = np.copy(self.ci)


    def norm_ci(self, *args):
        self.ci = self.ci / np.sqrt(np.dot(np.conj(self.ci), self.ci))

















