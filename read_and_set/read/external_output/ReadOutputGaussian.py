import numpy as np
import pandas as pd


### e se chiamassimo questa ReadExternalOutput e ci mettiamo tutto?


class ReadOutputGaussian():

    def read_ci0(self, name_file):
        ci0 = np.loadtxt(name_file)
        return ci0

    def read_en_ci0(self, name_file):
        en_ci0 = np.loadtxt(name_file, usecols=(3,))
        en_ci0 = np.insert(en_ci0, 0, 0)
        en_ci0 = en_ci0/27.2114
        return en_ci0

    def read_muT(self, name_file, n_en):
        muT_file = np.loadtxt(name_file, usecols=(4, 5, 6))
        muT = self.read_half_below_matrix(n_en, 3, muT_file)
        return muT

    def read_muT_gamess(self, name_file, n_en):
        file = np.loadtxt(name_file, usecols=(4, 5, 6))
        muT = self.read_half_above_matrix(n_en, 3, file)
        return muT

    def read_N_tessere_cavity(self, name_file): #read N tessere from a file which has N tessere on the first row
        n_tessere_cavity=pd.read_csv(name_file, nrows=1, header=None,engine='python')
        return n_tessere_cavity.iloc[0, 0]


    def read_V(self, name_file, n_en): #n_en is total number of states, nexcited+1
        print(n_en)
        n_tessere_cavity = self.read_N_tessere_cavity(name_file)
        VN = pd.read_csv(name_file, header=None, skiprows=2, nrows=n_tessere_cavity, sep=r"\s+", engine='python')
        VN = self.convert_fortran_dble(VN[2])
        index=list()
        index.append(n_tessere_cavity)
        n = int(n_en*(n_en+1)/2) # total number of elements in gaussian format for triangular mat: 00 01 ..0n 11 21 22 31...
        for i in range(n-2):
            index.append(index[-1] + n_tessere_cavity + 1) #identify headers. index[0]=n_tessere then always n_tessere+1
        index = np.asarray(index)
        V = pd.read_csv(name_file, header=None, skiprows=2, sep=r"\s+", engine='python')
        V = V.drop(V.index[index]) #dropping separation header lines
        V = self.convert_fortran_dble(V[0])
        V = V.reshape((n, n_tessere_cavity))
        V_ijn_el = self.read_half_below_matrix(n_en, n_tessere_cavity, V)
        V_tot = -np.array(V_ijn_el)
        for i in range(n_en):
            V_tot[i,i,:] = -V_tot[i,i,:] + VN
        return V_tot


    #read compare wavet gamess format
    def read_V_wavet(self, name_file, n_en): #n_en is total number of states, nexcited+1
        print(n_en)
        n_tessere_cavity = self.read_N_tessere_cavity(name_file)
        VN = pd.read_csv(name_file, header=None, skiprows=2, nrows=n_tessere_cavity, sep=r"\s+", engine='python')
        VN = self.convert_fortran_dble(VN[2])
        index=list()
        index.append(n_tessere_cavity)
        n = int(n_en*(n_en+1)/2) # total number of elements in gaussian format for triangular mat: 00 01 ..0n 11 21 22 31...
        for i in range(n-2):
            index.append(index[-1] + n_tessere_cavity + 1) #identify headers. index[0]=n_tessere then always n_tessere+1
        index = np.asarray(index)
        V = pd.read_csv(name_file, header=None, skiprows=2, sep=r"\s+", engine='python')
        V = V.drop(V.index[index]) #dropping separation header lines
        V = self.convert_fortran_dble(V[0])
        V = V.reshape((n, n_tessere_cavity))
        V_ijn_el = self.read_half_below_matrix(n_en, n_tessere_cavity, V)
        V_tot = np.array(V_ijn_el)
        for i in range(n_en):
            V_tot[i, i, :] = V_tot[i, i, :] + VN
        return V_tot


    def read_V_gamess_noVN_sottr(self, name_file, n_en):
        n_tessere_cavity =self.read_N_tessere_cavity(name_file)
        index = list()
        index.append(n_tessere_cavity)
        n = int(n_en*(n_en+1)/2)
        for i in range(n-2):
            index.append(index[-1] + n_tessere_cavity + 1)
        index = np.asarray(index)  #index allowing to identify separation header lines
        V = pd.read_csv(name_file, header=None, skiprows=2, sep=r"\s+", engine='python')
        V = V.drop(V.index[index]) #dropping separation header lines
        V = self.convert_fortran_dble(V[0])
        V = V.reshape((n, n_tessere_cavity))
        V_ijn_el = self.read_half_below_matrix(n_en, n_tessere_cavity, V)
        V_tot = np.array(V_ijn_el)
        return V_tot

    #funziona!!!
    def convertVgamess_to_VWaveT(self, namein, nameout, n_en):
        n_tessere_cavity =self.read_N_tessere_cavity(namein)
        V = pd.read_csv(namein, header=None, skiprows=2, sep=r"\s+", engine='python')
        n_0_pot = n_tessere_cavity * n_en + n_en-1
        V0 = np.array(V[:n_0_pot])
        f = open(nameout, 'w+')
        V_gamess_exc = self.read_V_gamess_noVN_sottr(namein, n_en)
        f.write(str(V_gamess_exc.shape[2]) + "\n")
        f.write("0   0  \n")
        f.close()
        f = open(nameout, 'ab')
        np.savetxt(f, V0, delimiter=' ',  fmt="%s")
        f.close()
        for i in range(1, n_en):
            for j in range(i + 1):
                if j != 0:
                    f = open(nameout, 'a')
                    f.write(str(i) + "  " + str(j) + "\n")
                    f.close()
                    f = open(nameout, 'ab')
                    np.savetxt(f, V_gamess_exc[i, j], delimiter=' ', header='', footer='', fmt='%1.8f')
                    f.close()



    def read_Q_matrix(self, Q_name_file):
        Q_tmp = pd.read_csv(Q_name_file, header=None, engine='python')
        Q_tmp = self.convert_fortran_dble(Q_tmp[0])
        Q_gamess = Q_tmp.reshape(int(np.sqrt(Q_tmp.size)), int(np.sqrt(Q_tmp.size)))
        return Q_gamess

    def read_cavity_tesserae(self, cav_file):
        n_tessere, n_skip = pd.read_csv(cav_file, nrows=1, header=None, sep="\s+").iloc[0,:]
        cavity = pd.read_csv(cav_file, header=None, skiprows=n_skip+1, sep=r"\s+")
        cavity = np.array(cavity)
        return cavity


    def convert_fortran_dble(self, pandas_data):
        converted = np.zeros_like(pandas_data)
        for i in range(converted.size):
            converted[i]=float(pandas_data.iloc[i].replace("D","E"))
        return converted

# read gaussian files with format: 00x 00y 00z, 01x 01y 01z...., 11x 11y 11z, 21x 21y 21z, 22x 22y 22z, 31x...
    def read_half_below_matrix(self, n_row, n_3D_col, half_matrix):
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



    def read_half_above_matrix(self, n_row, n_3D_col, half_matrix):
        matrix = np.zeros(shape=(n_row, n_row, n_3D_col))
        trans_0i = half_matrix[0:n_row]
        tmp_ij = half_matrix[n_row:]
        iu = np.triu_indices(n_row-1)
        trans_ij = np.zeros(shape=(n_row-1, n_row-1, n_3D_col))
        trans_ij[iu] = tmp_ij
        matrix[0] = trans_0i
        matrix[1:, 1:] = trans_ij
        matrix[0, :] = matrix[:, 0]
        for i in np.arange(1, n_row):
            matrix[:, i] = matrix[i, :]
        return matrix

