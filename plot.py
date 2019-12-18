import numpy as np
import matplotlib.pyplot as plt
from sys import exit

N = 100
L = 17.997
rho = N/(L**3)

A = 1.
alpha = 1.8

def g1(rho):
    eta = np.pi*rho/6
    return (1+0.5*eta)/( (1-eta)**2 )


def phi(r):
    u = -A*np.exp(-alpha*(r-1) )/r
    return np.exp(-u)
rr = np.linspace(1,10,1000) 
grr = phi(rr)

gr = np.loadtxt("gr.dat")
r = gr[:,0] + 1
g = gr[:,1]

gr = np.loadtxt("gr1.dat")
r1 = gr[:,0] + 1
g1 = gr[:,1]


gr = np.loadtxt("sim.dat")
rsim = gr[:,0]
gsim= gr[:,1]

#plt.plot(r1,g1, color='blue', label="square")
plt.plot(r,g, color='black', label = "shift")
#plt.plot(rsim,gsim, color='red', label = "data")
#plt.plot(rr,grr, color='red')


plt.title(rho)
plt.show()

