import numpy as np
from save.Save import Save
import pandas as pd
from read import auxiliary_functions as af
from pcm.PCM import PCM, FrozenSolventPCM, DinamicPCM

class Env():
    def __init__(self):
        self.env = None
        self.pcm = PCM()
        self.muLF = 0
        self.save = Save()

    def set_env(self, env):
        self.env = env

    def get_env(self):
        return self.env

    def get_q_t(self):
        #q_tmp= self.q_t.q_t[0]+self.q_t.q_t[1]
        #diff= 0 - q_tmp
        #print(diff)
        #delta_x_tessera = diff/136
        #q_tot = q_tmp + delta_x_tessera
        #return q_tot
        return self.pcm.q_t[0] + self.pcm.q_t[1]



    def get_q_t_lf(self):
        return self.pcm.qijn_lf

    def set_q_t(self, PCM_parameters, mol, field_t):
        if self.env == "vac":
            pass
        else:
            self.pcm.init_pcm(PCM_parameters, mol, field_t)
        return

    def propagate_q_t(self, mol, field_t):
        if self.env == "vac":
            pass
        else:
            self.pcm.propagate(mol, field_t)


    def propagate_q_t_bwd_oc(self, mol, field_t, ci):
        if self.env == "vac":
            pass
        else:
            self.pcm.propagate_bwd_oc(mol, field_t, ci)






class EnvPy(Env):
    def __init__(self):
        super().__init__()
        self.q_t = FrozenSolventPCM()

    def set_muLF(self, mol):
        self.muLF = -af.matrix_prod_tesserae_ijn_nn(self.q_t.qijn_lf, mol.Vijn)






class EnvTDPlas(Env):
    def __init__(self):
        super().__init__()
        self.q_t = DinamicPCM()








class PCM():
    def __init__(self):
        super().__init__()
        self.env = None
        self.Vij = None
        self.cavity = None
        self.qijn = None
        self.qijn_flip = None
        self.Vij_flip = None
        self.muT_local_field = None
        self.to_subtract_fromH = None
        self.save = Save()



    def set_environment(self, env):
        self.env = env

    def initialize_pcm(self, read_qijn, name, Vij, Q_gamess, cavity, Q_local_field_gamess, h):
        self.set_cavity(cavity)
        self.set_Vij(Vij)
        self.set_qijn(read_qijn, name, Q_gamess)
        self.calc_muT_local_field(Q_local_field_gamess, h)
        self.set_to_subtract_fromH()
        self.qijn_flip = af.flip_3D_py2f(self.qijn)
        self.Vij_flip = af.flip_3D_py2f(self.Vij)



    # <editor-fold desc="SetsGets">
    def set_Vij(self, Vij):
        self.Vij = Vij

    def set_cavity(self, cavity):
        self.cavity = cavity

    def set_qijn(self, read_qijn, name, Q_gamess):
        if(read_qijn == True):
            self.read_qijn(name)
        else:
            self.calc_qijn(name, Q_gamess)

    def set_to_subtract_fromH(self):
        to_subtract_fromH = np.zeros([self.qijn.shape[0], self.qijn.shape[0]])
        for n in range(self.qijn.shape[2]):
            to_subtract_fromH[:,:] += self.qijn[0, 0, n]*self.Vij[:,:,n]
        self.to_subtract_fromH = to_subtract_fromH
        return to_subtract_fromH

    def get_env(self):
        return self.env

    def get_Vij(self):
        return self.Vij

    def get_qijn(self):
        return self.qijn
# </editor-fold>


# <editor-fold desc="calculate">
    def calc_muT_local_field(self, Q_local_field_gamess, h):
        Q_local_field_tdplas = np.dot(Q_local_field_gamess, np.diag(self.cavity[:,3]))
        q_local = - np.dot(Q_local_field_tdplas, self.cavity[:,0:3])
        h.mod_muT_local_field(-af.matrix_prod_tesserae_ijn_nn(q_local, self.Vij))


    def calc_qijn(self, name, Q_gamess):
        Q_tdplas = np.dot(Q_gamess, np.diag( self.cavity[:,3]))
        self.qijn = af.matrix_prod_tesserae_ijn_nn(Q_tdplas, self.Vij)

    def calc_and_write_qijn(self, name, Q_gamess):
        Q_tdplas = np.dot(Q_gamess, np.diag( self.cavity[:,3]))
        self.qijn = np.zeros(self.Vij.shape)
        self.save.create_bkp_file(name)
        f = open(name, 'w')
        f.write(str(Q_tdplas.shape[0])+"\n")
        for n in range(Q_tdplas.shape[0]):
            tmp = 0
            for h in range(Q_tdplas.shape[0]):
                tmp += Q_tdplas[n,h]*self.Vij[:,:,h]
            f.write("Cavity tessera "+str(n+1)+"\n")
            np.savetxt(f, tmp)
            self.qijn[:,:,n] = tmp
        f.close()

# </editor-fold>

    def read_qijn(self,name_file):
        n_tessere= pd.read_csv(name_file, nrows=1, header=None).iloc[0,0]
        index = []
        index.append(self.Vij.shape[0])
        for i in range(n_tessere-2):
            index.append(index[-1] + self.Vij.shape[0] + 1)
        index = np.asarray(index[1:])
        q = pd.read_csv(name_file, header=None, skiprows=2,sep=r"\s+")
        q = q.drop(q.index[index])
        q = np.asarray(q)
        q = q.astype(float(q))
        q.reshape(self.Vij[2], self.Vij[0], self.Vij[1])
        self.qijn = q.transpose((1,2,0))













