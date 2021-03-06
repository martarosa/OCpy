from medium.ABCMedium import ABCMedium
from read_and_set.read import auxiliary_functions as af
import numpy as np
from medium import interface_tdplas as tdplas
from parameters.MediumParameters import MediumParameters

#When the medium is treated as dinamicaly changing during the propagation, the charges on the solvent cavity
#or nanoparticle tessere are propagated in TDPlas
#OCpy gives to TDPlas the coefficients and the external field at each time step and receives back
#the propagated charges

class DinamicMedium(ABCMedium):
    def __init__(self):
        super().__init__()
        self.par = MediumParameters()
        #self.internal_field = [] local field NP



    def init_medium(self, medium_input, mol, field_object):
        self.par.medium = medium_input.medium
        self.par.cavity = medium_input.cavity
        Vijn_mol_fortran_flip = af.flip_3D_py2f(mol.par.Vijn)
        self.par.q_t = np.zeros(mol.par.Vijn.shape[2])
        tdplas.interface_tdplas.init_global_tdplas(field_object.time_axis[1],#in pratica gli sto passando il dt
                                  Vijn_mol_fortran_flip,
                                  Vijn_mol_fortran_flip.shape[0],
                                  Vijn_mol_fortran_flip.shape[1])
        tdplas.interface_tdplas.call_init_charges(mol.wf.ci_prev[0], field_object.f_xyz[0])



    def reset_medium(self, mol, field_object):
        self.par.q_t = np.zeros(mol.par.Vijn.shape[2])
        tdplas.interface_tdplas.call_init_charges(mol.wf.ci_prev[0], field_object.f_xyz[0])
        #self.get_q_t() debugging di TDPLAS

    def propagate_charges(self, mol, field_dt_vector):
        af.exit_error("ERROR: Propagation without TDPLAS is not possible for the nanoparticle")



    def propagate_charges_fortran(self, mol, field_dt_vector):
        tdplas.interface_tdplas.call_prop_charges(mol.wf.ci_prev[0],
                                                  field_dt_vector,
                                                  self.par.q_t)

    def get_q_t(self):
        return self.par.q_t


    #prit local field for np case
    #def calc_local_field_t(self):
    #    mol_cc = np.array([24,0.0,0.0])
    #    cavity = self.par.cavity
    #    diff = np.zeros(shape=(3))
    #    internal_field_t = np.zeros(shape=(3))
    #    for i in range(cavity.shape[0]):
    #        diff[0] = mol_cc[0] - cavity[i, 0]
    #        diff[1] = mol_cc[1] - cavity[i, 1]
    #        diff[2] = mol_cc[2] - cavity[i, 2]
    #        dist = np.sqrt(np.dot(diff,diff))
    #        internal_field_t = internal_field_t + self.par.q_t[i]*diff/(dist*dist*dist)
    #    self.internal_field.append(internal_field_t)


