import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

beta = 1

def Psi(x, alpha):
    A = ((8*alpha**3)/(beta*np.pi**3))**(1/4)
    return A*np.exp(-alpha*x*x)

x = np.linspace(0, 6, 1000)
label_size = {"size":"14"}
sns.set()

plt.plot(x, Psi(x, 0.5), '-',  label=r'$\alpha=0.5$')
plt.plot(x, Psi(x, 0.35), '--', label=r'$\alpha=0.35$')
plt.plot(x, Psi(x, 0.2), '-.', label=r'$\alpha=0.2$')
plt.xlabel(r"r [a$_{ho}$]", **label_size)
plt.ylabel("Wave function", **label_size)
plt.legend()
#plt.savefig('../images/wavefunctions.png')
plt.show()
