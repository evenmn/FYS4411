import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os

# Numerical values
dirpath = os.getcwd()
data = np.loadtxt("../data/ob_density.dat")
#data1 = np.loadtxt("../data/ob_density_a_0.dat")
label_size = {"size":"14"}

r = np.linspace(0, 5, len(data))

N = 100     # Number of particles

def f(x):
    return (1.8*N/(8*np.pi**(3/2)))*np.exp(-(x**2))

sns.set()
#plt.plot(r, data1, linewidth=2.0, label="a=0.0 Numerical")
plt.plot(r, f(r), '--r', linewidth=2.0, label="a=0.0 Analytical")
plt.plot(r, data, linewidth=2.0, label="a=0.0043")
#plt.title("One Body Density", **label_size)
plt.xlabel("r", **label_size)
plt.ylabel(r"$\rho$", **label_size)
plt.legend(loc="best")
#plt.savefig("%s/images/ob.png" % dirpath)
plt.show()
