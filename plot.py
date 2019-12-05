import numpy as np
import matplotlib.pyplot as plt
from sys import exit

gr = np.loadtxt("gr.dat")
r = gr[:,0]
g = gr[:,1]

plt.plot(r,g)
plt.show()

