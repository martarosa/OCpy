import numpy as np

from abc import ABCMeta, abstractmethod

from SystemObj import DiscreteTimePar
from field.Field import Func_tMatrix
from parameters.OCIteratorParameters import OCIteratorParameters
from ObjectiveFunction import ObjectiveFunction


class ABCOCIterator(metaclass=ABCMeta):
    def __init__(self):
        self.obj_fun = ObjectiveFunction()
        self.par = OCIteratorParameters()
        self.discrete_t_par = DiscreteTimePar()

        self.field_psi_matrix = Func_tMatrix()
        self.psi_coeff_t_matrix = Func_tMatrix()

        self.dict_out = {}


    @abstractmethod
    def iterate(self, current_iteration):
        pass

    @abstractmethod
    def check_convergence(self):
        pass

    @abstractmethod
    def calc_J(self):
        pass

    @abstractmethod
    def init(self, molecule, starting_field, medium, alpha_t, oc_input, oc_conf, prop_conf):
        pass

    @abstractmethod
    def init_output_dictionary(self):
        pass


    @abstractmethod
    def get_restart(self):
        pass


    def field_J_integral(self):
        ax_square = self.field_psi_matrix.f_xyz.ndim - 1
        ax_integral= self.field_psi_matrix.f_xyz.ndim - 2
        f_square = np.sum(self.field_psi_matrix.f_xyz * self.field_psi_matrix.f_xyz, axis=ax_square)
        f_integral = np.sum(f_square, axis=ax_integral)
        out_field = f_integral*self.discrete_t_par.dt
        return out_field


    def alpha_field_J_integral(self):
        ax_square= self.field_psi_matrix.f_xyz.ndim - 1
        ax_integral= self.field_psi_matrix.f_xyz.ndim - 2
        f_square = np.sum(self.field_psi_matrix.f_xyz * self.field_psi_matrix.f_xyz, axis=ax_square) * self.par.alpha_t
        f_integral = np.sum(f_square, axis=ax_integral)
        out_integral = f_integral*self.discrete_t_par.dt
        return out_integral