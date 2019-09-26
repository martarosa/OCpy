import numpy as np
import time
import matmult as mm 


def double_summation(ket, bra, M):
    double_sum = np.dot(bra, np.sum(M*ket[np.newaxis, :, np.newaxis], axis=1))
    return double_sum

def single_summation_tessere(ket, M):
    single_sum = np.sum(M*ket[np.newaxis, np.newaxis, :], axis=2)
    return single_sum


M=np.array([[[1.1+1j*1.1,2.1+1j*2.1],[3.1+1j*3.1,4.1+1j*4.1]],[[5.1+1j*5.1,6.1+1j*6.1],[7.1+1j*7.1,8.1+1j*8.1]],[[9.1+1j*9.1,10.1+1j*10.1],[11.1+1j*11.1,12.1+1j*12.1]]])

Mpy=np.array([[[1.1+1j*1.1,5.1+1j*5.1,9.1+1j*9.1],[2.1+1j*2.1,6.1+1j*6.1,10.1+1j*10.1]],[[3.1+1j*3.1,7.1+1j*7.1,11.1+1j*11.1],[4.1+1j*4.1,8.1+1j*8.1,12.1+1j*12.1]]])

ket=np.array([1,2])
bra=np.array([-1,-2])
M1=np.array([[3.1+1j*3.1,4.1+1j*4.1],[4.1+1j*4.1,5.1+1j*5.1]])
M2=np.array([[-3.1+1j*3.1,4.1+1j*4.1],[-4.1+1j*4.1,5.1+1j*5.1]])
Mr=np.array([[3.1,4.1],[4.1,5.1]])
#M=np.asfortranarray(M)



#prod = mm.doublesum(ket,bra,M,)
#prod2 = mm.single_summation_tessere(prod,M)

print(Mr.dtype)

prod=mm.dot(np.asfortranarray(ket, dtype=np.complex64),np.asfortranarray(M2, dtype=(np.complex64)))


#prodpy=double_summation(ket,bra,Mpy)
#prodpy2=single_summation_tessere(prodpy,Mpy)
#print(prod2[1])
#print(prodpy2[1])
