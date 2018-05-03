import numpy as np
import matplotlib.pyplot as plt

asymptote = 0.5

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

'''
data1 = np.loadtxt("../data/BF_1P_N_1_MC_1000000_eta_001.txt")
data2 = np.loadtxt("../data/BF_1P_N_2_MC_1000000_eta_001.txt")
data3 = np.loadtxt("../data/BF_1P_N_3_MC_1000000_eta_001.txt")
data5 = np.loadtxt("../data/BF_1P_N_5_MC_1000000_eta_001.txt")
data6 = np.loadtxt("../data/BF_1P_N_10_MC_1000000_eta_001.txt")

x1 = np.linspace(0, len(data1) - 1, len(data1))
x2 = np.linspace(0, len(data2) - 1, len(data2))
x3 = np.linspace(0, len(data3) - 1, len(data3))
x5 = np.linspace(0, len(data5) - 1, len(data5))
x6 = np.linspace(0, len(data6) - 1, len(data6))

plt.plot(x1, data1, label="N=1")
plt.plot(x2, data2, label="N=2")
plt.plot(x3, data3, label="N=3")
plt.plot(x5, data5, label="N=5")
plt.plot(x6, data6, label="N=6")
plt.axhline(asymptote, linestyle='--', color='r', label="Exact")
plt.xlabel("Iteration")
plt.ylabel("Energy")
plt.grid()
plt.legend()
#plt.axis([-10, 210, 1.5, 6])
#plt.savefig('figure.png')
plt.show()
'''
