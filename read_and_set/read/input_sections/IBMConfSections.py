from read_and_set.read.input_sections.ABCSection import ABCSection


class SectionIBMParameters(ABCSection):
    def __init__(self):
        super().__init__()
        self.section = "IBM_INTERFACE"
        self.section_default_dictionary = {'provider' : 'qasm_simulator',
                                           'device' : 'none',
                                           'noise' : 'False',
                                           'mapping' : 'linear',
                                           'shots' : '5000'}
        self.section_dictionary = {}
        self.allowed_val = [['provider', ['qasm_simulator', 'statevector_simulator']],
                            ['device', ['ibmq_santiago', 'ibmq_16_melbourne','ibmqx2']],
                            ['mapping',['linear', 'binary']]
                            ['noise', ['True', 'False']]]
        self.case_unsensitive_keys = ['provider', 'device']