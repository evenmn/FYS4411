import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os

#Analytical values
def E_L(alpha):
    return alpha

# Numerical values
dirpath = os.getcwd()
data = np.loadtxt("%s/data/local_energy.dat" % dirpath)
label_size = {"size":"14"}

sns.set()
plt.plot(data[:,0], data[:,1], label="Numerical")
plt.title("Local energy", **label_size)
plt.xlabel(r"$\alpha$", **label_size)
plt.ylabel(r"$E_L(\alpha)$", **label_size)
plt.legend(loc="best")
plt.savefig("%s/images/energy.png" % dirpath)
plt.show()

plt.plot(data[:,0], data[:,2], label="Numerical")
plt.title("Variance", **label_size)
plt.xlabel(r"$\alpha$", **label_size)
plt.ylabel(r"$\sigma$", **label_size)
plt.legend(loc="best")
plt.savefig("%s/images/variance.png" % dirpath)
