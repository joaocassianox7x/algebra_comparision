import numpy as np
import matplotlib.pyplot as plt
import time

z=np.loadtxt("tempos.dat")

ORD = 25
ORD_N = ORD*100
z1 = np.zeros(ORD)
t = np.zeros(ORD)
t1 = np.zeros(ORD)

for i in range(1,ORD+1):

    t0=time.time()
    np.linalg.eig(np.random.rand(i*100,i*100))
    t[i-1] = time.time()-t0
    
    t0=time.time()
    np.linalg.eigh(np.random.rand(i*100,i*100))
    t1[i-1] = time.time()-t0
    


plt.figure()
x=np.arange(ORD)*100+100


plt.plot(x,t,'r-1')
plt.plot(x,t1,'k->')

plt.plot(x,z[:,0],'y-o')

plt.plot(x,z[:,1],'b--')

plt.plot(x,z[:,2],'g-x')
plt.plot(x,z[:,3],'c-')


plt.legend(["NUMPY EIG","NUMPY EIGH","LAPACK DSTEVD", "LAPACK DSYEV","LAPACK DGEEV","ARMADILLO"])
plt.xlabel("Ordem (n)")
plt.ylabel("Tempo (s)")
plt.xlim(100,ORD_N)

plt.savefig("performance.png")