import numpy as np
import matplotlib.pyplot as plt
from sys import exit

rho = 0.6
def g1(rho):
    eta = np.pi*rho/6
    return (1+0.5*eta)/((1-eta)**2)
def h1(rho):
    return -1+g1(rho)


data = np.loadtxt("rho06_pc.dat")
plt.plot(data[:,0], data[:,1])
plt.scatter(1,g1(rho))
plt.axhline(1., color='black')
plt.show()

