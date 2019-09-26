import numpy as np


class ReadMethods:

    def convert_fortran_dble(self, pandas_data):
        converted = np.zeros_like(pandas_data)
        for i in range(converted.size):
            converted[i]=float(pandas_data.iloc[i].replace("D","E"))
        return converted

# read gaussian files with format: 00x 00y 00z, 01x 01y 01z...., 11x 11y 11z, 21x 21y 21z, 22x 22y 22z, 31x...
    def read_half_above_matrix_gaussian(self, n_row, n_3D_col, half_matrix):
        matrix = np.zeros(shape=(n_row, n_row, n_3D_col))
        trans_0i = half_matrix[0:n_row]
        tmp_ij = half_matrix[n_row:]
        il = np.tril_indices(n_row-1)
        trans_ij = np.zeros(shape=(n_row-1, n_row-1, n_3D_col))
        trans_ij[il] = tmp_ij
        matrix[0] = trans_0i
        matrix[1:, 1:] = trans_ij
        matrix[:, 0] = matrix[0, :]
        for i in np.arange(1, n_row):
            matrix[i, :] = matrix[:, i]
        return matrix


    def read_half_below_matrix_gaussian(self, n_row, n_3D_col, half_matrix):
        matrix = np.zeros(shape=(n_row, n_row, n_3D_col))
        trans_0i = half_matrix[0:n_row]
        tmp_ij = half_matrix[n_row:]
        iu = np.triu_indices(n_row-1)
        trans_ij = np.zeros(shape=(n_row-1, n_row-1, n_3D_col))
        trans_ij[iu] = tmp_ij
        matrix[0] = trans_0i
        matrix[1:, 1:] = trans_ij
        matrix[:, 0] = matrix[0, :]
        for i in np.arange(1, n_row):
            matrix[:, i] = matrix[i, :]
        return matrix