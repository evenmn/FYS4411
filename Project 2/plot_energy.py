import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("energy.txt")

x = np.linspace(0, len(data) - 1, len(data))

plt.plot(x, data)
plt.axhline(2.0, linestyle='--', color='r')
plt.grid()
plt.show()
