import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os

# Numerical values
dirpath = os.getcwd()
data0 = np.loadtxt("../data/ob_density_a_0_N_100.dat")
data1 = np.loadtxt("../data/ob_density_a_00043_N_100.dat")
data2 = np.loadtxt("../data/ob_density_a_0043_N_100.dat")
#data3 = np.loadtxt("../data/ob_density_a_043.dat")
label_size = {"size":"14"}

x = np.linspace(0, 3, len(data0))
r = np.sqrt(3*x**2)

N = 10     # Number of particles

def f(x, y, z, alpha, beta):
    return (2*alpha/np.pi)**(3/2)*(beta)**(1/2)*np.exp(-2*alpha*(x**2 + y**2 + beta*z**2))

sns.set()
plt.plot(r, data0, '-', linewidth=1.0, label="a=0.0")
#plt.plot(r, f(x, x, x, alpha=0.5, beta=np.sqrt(8)), '--r', linewidth=2.0, label="a=0.0 Analytical")
plt.plot(r, data1, '-',linewidth=1.0, label="a=0.0043")
plt.plot(r, data2, '-', linewidth=1.0, label="a=0.043")
#plt.plot(r, data3, '-', linewidth=1.0, label="a=0.43")
#plt.title("One Body Density", **label_size)
plt.xlabel("r$[a_{ho}]$", **label_size)
plt.ylabel(r"$\rho$", **label_size)
plt.legend(loc="best")
#plt.savefig("%s/images/ob.png" % dirpath)
plt.show()
