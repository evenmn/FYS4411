import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os

# Numerical values
dirpath = os.getcwd()
data = np.loadtxt("data/ob_density.dat")
label_size = {"size":"14"}

r = np.linspace(0, 5, len(data))

sns.set()
plt.plot(r, data, label="Numerical")
plt.title("One Body Densities, $a=1.0$", **label_size)
plt.xlabel("r", **label_size)
plt.ylabel(r"$\rho$", **label_size)
#plt.legend(loc="best")
plt.savefig("%s/images/ob_a_1.png" % dirpath)
plt.show()

