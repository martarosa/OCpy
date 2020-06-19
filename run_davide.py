import argparse
import sys
sys.path.append('/scratch/martarosa/python/26-02-2020/OCpy/')
import SystemManager as ini
import time


folder = "./"
namefile = "input.dat"


from qiskit import QuantumRegister, ClassicalRegister, QuantumCircuit

from qiskit import Aer, execute

qreg = QuantumRegister(2, 'q')

creg = ClassicalRegister(2, 'c')

circuit = QuantumCircuit(qreg, creg)

circuit.x(qreg[0])
circuit.x(qreg[1])

provider = Aer.get_backend("statevector_simulator")

result = execute(circuit, provider).result()

statevector = result.get_statevector(circuit)

print(statevector)


OC_system=ini.SystemManager()
OC_system.init_system(folder, namefile)



#start=time.time()
OC_system.oc.iterate()
#end=time.time()
#print("serial: " + str(end-start))

