import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os

#Analytical values
def EGP(alpha, beta, N):
    return N*(1 + beta/2)*(1/(4*alpha) + alpha)

# Numerical values
dirpath = os.getcwd()
data = np.loadtxt("%s/../data/local_energy.dat" % dirpath)
label_size = {"size":"14"}

alpha = np.linspace(0.1, 0.9, 1000)

sns.set()
plt.plot(alpha, EGP(alpha, 1, 10), label="GP")
plt.plot(data[:,0], data[:,1], label="VMC")
#plt.title("Local energy", **label_size)
plt.xlabel(r"$\alpha$", **label_size)
plt.ylabel(r"$E_L(\alpha)$", **label_size)
plt.legend(loc="best")
plt.savefig("%s/../images/energy.png" % dirpath)
plt.show()

plt.plot(data[:,0], data[:,2], label="Numerical")
#plt.title("Variance", **label_size)
plt.xlabel(r"$\alpha$", **label_size)
plt.ylabel(r"$\sigma$", **label_size)
#plt.legend(loc="best")
plt.savefig("%s/../images/variance.png" % dirpath)
plt.show()
