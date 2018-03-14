import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os

# Numerical values
dirpath = os.getcwd()
data = np.loadtxt("data/ob_density.dat")
label_size = {"size":"14"}

r = np.linspace(0, 5, len(data))

N = 10     # Number of particles

def f(x):
    return (np.sqrt(np.pi)/2)**(N-1)*np.exp(-(x**2))

sns.set()
plt.plot(r, data, label="Numerical")
plt.plot(r, f(r), label="Analytical")
#plt.title("One Body Densities, $a=1.0$", **label_size)
plt.xlabel("r", **label_size)
plt.ylabel(r"$\rho$", **label_size)
plt.legend(loc="best")
#plt.savefig("%s/images/ob_a_1.png" % dirpath)
plt.show()

