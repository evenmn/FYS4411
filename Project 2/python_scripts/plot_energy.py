import numpy as np
import matplotlib.pyplot as plt

asymptote = 2.0

data = np.loadtxt("../data/energy.txt")
x = np.linspace(0, len(data) - 1, len(data))

plt.plot(x, data)
plt.axhline(asymptote, linestyle='--', color='r')
plt.grid()
#plt.savefig('figure.png')
plt.show()
