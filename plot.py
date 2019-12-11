import numpy as np
import matplotlib.pyplot as plt
from sys import exit

N = 1000
L = 20.
rho = N/(L**3)

def g1(rho):
    eta = np.pi*rho/6
    return (1+0.5*eta)/( (1-eta)**2 )


gr = np.loadtxt("gr.dat")
r = gr[:,0]
g = gr[:,1]

plt.plot(r,g)
plt.scatter(1,g1(rho) )
plt.title(rho)
plt.show()

