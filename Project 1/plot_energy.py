import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os

#Analytical values
def E_L(alpha):
    return alpha

# Numerical values
dirpath = os.getcwd()
print(dirpath)
data = np.loadtxt("%s/data/local_energy.dat" % dirpath)
print(data)

sns.set()
plt.plot(data[:,0], data[:,1], label="Numerical")
plt.title("Local energy")
plt.xlabel(r"$\alpha$")
plt.ylabel(r"$E_L(\alpha)$")
plt.legend(loc="best")
plt.savefig("%s/images/energy.png" % dirpath)
plt.show()

plt.plot(data[:,0], data[:,2], label="Numerical")
plt.title("Variance")
plt.xlabel(r"$\alpha$")
plt.ylabel(r"$\sigma$")
plt.legend(loc="best")
plt.savefig("%s/images/variance.png" % dirpath)
plt.show()
