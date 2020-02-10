from pcm.ABCPCM import ABCPCM

# if the solvent is dinamic qijn and qijn_lf depend on time = qijn_t, qijn_lf_t
#q_t dependence on time comes from <psi(t)|qijn_t|psi(t)>
#q_t_lf dependence comes from qijn_lf_t * field_t
#at each timestep qijn_t and qijn_lf_t are given by tdplas

class DinamicPCM(ABCPCM):
    def __init__(self):
        super().__init__()


    def init_pcm(self, PCM_input, mol, field_t):
        mana = 1
        #input_V_tdplas = af.double_summation(mol.wf.ci_prev[0], np.conj(mol.wf.ci_prev[0]), mol.Vijn)
        #input_V_lf_tdplas = - np.dot(field_t, self.cavity[:,0:3])
        #tdplas.init_medium(
        #                  np.asfortranarray(input_V_tdplas, dtype=np.float64),
        #                  np.asfortranarray(input_V_lf_tdplas, dtype=np.float64))


    def propagate(self, i, mol, field_t):
        mana  = 1
        #input_V_tdplas = af.double_summation(mol.wf.ci_prev[0], np.conj(mol.wf.ci_prev[0]), mol.Vijn)
        #input_V_lf_tdplas = - np.dot(field_t, self.cavity[:,0:3])
        #tdplas.propagate_mdm(i,
        #                  np.asfortranarray(input_V_tdplas, dtype=np.float64),
        #                  np.asfortranarray(input_V_lf_tdplas, dtype=np.float64))


    def get_q_t(self):
        # q_tmp= self.q_t.q_t[0]+self.q_t.q_t[1]
        # diff= 0 - q_tmp
        # print(diff)
        # delta_x_tessera = diff/self.q00n.shape[-1]
        # q_tot = q_tmp + delta_x_tessera
        # return q_tot
        return self.q_t[0] + self.q_t[1]