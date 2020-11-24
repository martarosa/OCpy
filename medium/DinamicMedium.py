from medium.ABCMedium import ABCMedium
from read_and_set.read import auxiliary_functions as af
import numpy as np
from medium import interface_tdplas as tdplas
from parameters.MediumParameters import MediumParameters
# if the solvent is dinamic qijn and qijn_lf depend on time = qijn_t, qijn_lf_t
#q_t dependence on time comes from <psi(t)|qijn_t|psi(t)>
#q_t_lf dependence comes from qijn_lf_t * field_dt_vector
#at each timestep qijn_t and qijn_lf_t are given by tdplas

class DinamicMedium(ABCMedium):
    def __init__(self):
        super().__init__()
        self.par = MediumParameters()



    def init_medium(self, medium_input, mol, field_object):

        self.par.medium = medium_input.medium
        Vijn_mol_fortran_flip = af.flip_3D_py2f(mol.par.Vijn)
        self.par.q_t = np.zeros(mol.par.Vijn.shape[2])
        tdplas.interface_tdplas.init_global_tdplas(field_object.time_axis[1],
                                  Vijn_mol_fortran_flip,
                                  Vijn_mol_fortran_flip.shape[0],
                                  Vijn_mol_fortran_flip.shape[1])
        tdplas.interface_tdplas.call_init_charges(mol.wf.ci_prev[0], field_object.f_xyz[0])



    def reset_medium(self, mol, field_object):
        self.par.q_t = np.zeros(mol.par.Vijn.shape[2])
        tdplas.interface_tdplas.call_init_charges(mol.wf.ci_prev[0], field_object.f_xyz[0])



    def propagate(self, mol, field_dt_vector):
        af.exit_error("ERROR: Propagation without TDPLAS is not possible for the nanoparticle")



    def propagate_fortran(self, mol, field_dt_vector):
        tdplas.interface_tdplas.call_prop_charges(mol.wf.ci_prev[0],
                                                  field_dt_vector,
                                                  self.par.q_t)




    def get_q_t(self):
        return self.par.q_t


