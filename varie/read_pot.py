import numpy as np
from read import ReadOutputGaussian as rgauss

folder="/home/mana/Desktop/python/optimal_control/quinolone/vac/"
namefile="input.dat"

nameV=folder+"ci_pot_gamess.inp"
outV=folder+"ci_pot.inp"
n_en=16

outG=rgauss.ReadOutputGaussian()
V=outG.read_V_gamess(nameV,n_en)

outG=rgauss.ReadOutputGaussian()
V=outG.read_V_gamess_noVN_sottr(nameV,n_en)


f=open(outV, 'w+')
f.write(str(V.shape[2])+"\n")
f.close()

for i in range(1,n_en):
    for j in range(i+1):
        if(j==0):
            print(j)
        else:
            f=open(outV,'a')
            f.write(str(i)+"  "+str(j)+"\n")
            f.close()
            f=open(outV,'ab')
            np.savetxt(f,V[i,j],delimiter = ' ', header = '', footer = '', fmt='%1.8f')
            f.close()


n05_diag=0
n05_no_diag=0

qij=np.sum(system.env.qijn, axis=2)
n_stati=qij.shape[0]


for i in range(n_stati):
    for j in range(n_stati):
        if(i==j):
            if(qij[i,j]>0.5) or (qij[i,j]<-0.5):
                n05_diag+=1
                print(str(i) +"  "+ str(j) +"  " +str(qij[i,j]))



for i in range(n_stati):
    for j in range(n_stati):
        if(i==j):
            print(i)
        else:
            if (qij[i, j] > 0.5) or (qij[i, j] < -0.5):
                print(str(i) + "  " + str(j) + "  " + str(qij[i, j]))
                n05_no_diag+=1
