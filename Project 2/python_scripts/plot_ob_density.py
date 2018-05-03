import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


# Numerical values
data1 = np.loadtxt("../data/OB_wo_interaction.dat")
data2 = np.loadtxt("../data/OB_w_interaction.dat")
label_size = {"size":"14"}

r1 = np.linspace(0, 3, len(data1))
r2 = np.linspace(0, 3, len(data2))

sns.set()
plt.plot(r1, data1, '-', linewidth=1.0, label="WO")
plt.plot(r2, data2, '-', linewidth=1.0, label="W")
plt.xlabel("r$[a_{ho}]$", **label_size)
plt.ylabel(r"$\rho$", **label_size)
plt.legend(loc="best")
plt.savefig("../plots/ob.png")
plt.show()
