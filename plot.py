import numpy as np
import matplotlib.pyplot as plt
from sys import exit

N = 1000
L = 20.
rho = N/(L**3)

def g1(rho):
    eta = np.pi*rho/6
    return (1+0.5*eta)/( (1-eta)**2 )

sim = np.loadtxt("sim_data_rho0.8_1.8_1.1111.dat")
simx = sim[:,0]
simg = sim[:,1]

gr = np.loadtxt("gr.dat")
r = gr[:,0]
g = gr[:,1]

gr1 = np.loadtxt("gr_rc3.dat")
r1 = gr[:,0]
g1 = gr[:,1]

plt.plot(r,g, color='black')
#plt.plot(r1,g1+0.1, color='blue')
plt.plot(simx,simg, color='red')
#plt.scatter(1,g1(rho) )
plt.title(rho)
plt.show()

