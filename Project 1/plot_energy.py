import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

data = np.loadtxt("local_energy.dat")

sns.set()
plt.plot(data[:,0], data[:,1], label="Numerical")
plt.title("Local energy")
plt.xlabel(r"$\alpha$")
plt.ylabel(r"$E_L(\alpha)$")
plt.legend(loc="best")
plt.show()
