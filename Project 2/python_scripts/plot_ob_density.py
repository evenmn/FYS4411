import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os

# Numerical values
dirpath = os.getcwd()
data = np.loadtxt("../data/ob_density.dat")
label_size = {"size":"14"}

r = np.linspace(0, 3, len(data))

sns.set()
plt.plot(r, data, '-', linewidth=1.0, label="a=0.0")
plt.xlabel("r$[a_{ho}]$", **label_size)
plt.ylabel(r"$\rho$", **label_size)
plt.legend(loc="best")
plt.savefig("../plots/ob.png")
plt.show()
