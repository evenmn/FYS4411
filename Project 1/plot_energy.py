import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

#Analytical values
def E_L(alpha):
    return alpha

# Numerical values
data = np.loadtxt("local_energy.dat")

sns.set()
plt.plot(data[:,0], data[:,1], label="Numerical")
plt.title("Local energy")
plt.xlabel(r"$\alpha$")
plt.ylabel(r"$E_L(\alpha)$")
plt.legend(loc="best")
plt.show()

plt.plot(data[:,0], data[:,2], label="Numerical")
plt.title("Variance")
plt.xlabel(r"$\alpha$")
plt.ylabel(r"$\sigma$")
plt.legend(loc="best")
plt.show()
