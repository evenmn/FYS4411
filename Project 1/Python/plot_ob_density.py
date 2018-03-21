import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os

# Numerical values
dirpath = os.getcwd()
data0 = np.loadtxt("../data/ob_density_a_0_N_10.dat")
data1 = np.loadtxt("../data/ob_density_a_00043_N_10.dat")
data2 = np.loadtxt("../data/ob_density_a_0043_N_10.dat")
label_size = {"size":"14"}

r = np.linspace(0, 3, len(data0))

N = 10     # Number of particles

def f(x):
    return np.pi**(-3/2)*8**(1/4)*np.exp(-x**2)

sns.set()
plt.plot(r, data0, linewidth=1.0, label="a=0.0")
#plt.plot(r, f(r), '--r', linewidth=2.0, label="a=0.0 Analytical")
plt.plot(r, data1, linewidth=1.0, label="a=0.0043")
plt.plot(r, data2, linewidth=1.0, label="a=0.043")
#plt.title("One Body Density", **label_size)
plt.xlabel("r$[a_{ho}]$", **label_size)
plt.ylabel(r"$\rho$", **label_size)
plt.legend(loc="best")
#plt.savefig("%s/images/ob.png" % dirpath)
plt.show()
