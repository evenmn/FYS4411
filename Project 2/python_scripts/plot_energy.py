import numpy as np
import matplotlib.pyplot as plt

asymptote = 3.0

data = np.loadtxt("../data/energy.txt")
x = np.linspace(0, len(data) - 1, len(data))

plt.plot(x, data, label="Calculated")
plt.axhline(asymptote, linestyle='--', color='r', label="Exact")
plt.xlabel("Iteration")
plt.ylabel("Energy")
plt.grid()
plt.legend()
#plt.axis([-10, 210, 1.5, 6])
#plt.savefig('figure.png')
plt.show()
