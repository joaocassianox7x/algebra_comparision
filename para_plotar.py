import numpy as np
import matplotlib.pyplot as plt

plt.figure()

z=np.loadtxt("output_potentials/vecx2.dat")
z1=np.loadtxt("output_potentials/valx2.dat")
x=np.linspace(0,1,1000)




for i in range(5):
    plt.plot(x,z1[i]+40*np.abs(z[:,i])**2/(x[1]-x[0]))
    plt.plot([0,1],[z1[i],z1[i]],"k")


plt.title(r"$V_0 = 2500x^2$")
plt.xlim(0,1)
plt.ylabel(r"$|\psi ^2|$")

plt.tight_layout()
plt.savefig("pot_x2.png")


plt.figure()

z=np.loadtxt("output_potentials/vec0.dat")
z1=np.loadtxt("output_potentials/val0.dat")
x=np.linspace(0,1,1000)




for i in range(5):
    plt.plot(x,z1[i]+40*np.abs(z[:,i])**2/(x[1]-x[0]))
    plt.plot([0,1],[z1[i],z1[i]],"k")


plt.title(r"$V_0 = 0$")
plt.xlim(0,1)
plt.ylabel(r"$|\psi ^2|$")

plt.tight_layout()
plt.savefig("pot_0.png")






plt.figure()
z=np.loadtxt("output_potentials/vecsin.dat")
z1=np.loadtxt("output_potentials/valsin.dat")
x=np.linspace(0,1,1000)




for i in range(5):
    plt.plot(x,z1[i]+40*np.abs(z[:,i])**2/(x[1]-x[0]))
    plt.plot([0,1],[z1[i],z1[i]],"k")


plt.title(r"$V = -1000sin(3 \pi x)$")
plt.xlim(0,1)
plt.ylabel(r"$|\psi ^2|$")

plt.tight_layout()
plt.savefig("pot_sin.png")








