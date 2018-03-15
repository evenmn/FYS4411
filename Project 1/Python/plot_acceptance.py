import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os

dirpath = os.getcwd()
dt = np.linspace(0.1, 1.9, 19)

# --- Importance Sampling ---
acceptance = np.array([0.999976, 0.999873, 0.999666, 0.999396,
                       0.998925, 0.998139, 0.996801, 0.994443,
                       0.990586, 0.985173, 0.977523, 0.967763,
                       0.955846, 0.941111, 0.923685, 0.903114,
                       0.879249, 0.852254, 0.821826])

sns.set()
label_size = {"size":"14"}

plt.plot(dt, acceptance, label="Numerical")
#plt.title("Acceptance rate Importance sampling", **label_size)
plt.xlabel("Timestep $\delta t$", **label_size)
plt.ylabel("Acceptance rate", **label_size)
#plt.legend(loc="best")
plt.savefig("%s/../images/acceptance_IS.png" % dirpath)
plt.show()

# --- Standard Metropolis ---
acceptance = np.array([0.971756, 0.943512, 0.915625, 0.888518,
                       0.859966, 0.833398, 0.806310, 0.779861,
                       0.754179, 0.729628, 0.704473, 0.680800,
                       0.656721, 0.634565, 0.613581, 0.591068,
                       0.569704, 0.550805, 0.532955])

sns.set()
label_size = {"size":"14"}

plt.plot(dt, acceptance, label="Numerical")
#plt.title("Acceptance rate brute force Metropolis", **label_size)
plt.xlabel("Step length $r$", **label_size)
plt.ylabel("Acceptance rate", **label_size)
#plt.legend(loc="best")
plt.savefig("%s/../images/acceptance_BF.png" % dirpath)
plt.show()
