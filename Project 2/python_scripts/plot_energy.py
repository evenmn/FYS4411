import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes, mark_inset

asymptote = 0.5
'''
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
data1 = np.loadtxt("../data/BF_1P_MC_1000000_eta_001_iter_500.txt")
data2 = np.loadtxt("../data/H_1P_MC_1000000_eta_001_iter_500_dt_01.txt")
data3 = np.loadtxt("../data/G_1P_N_1_MC_2pow20_eta_001_factor_1.txt")


x = np.linspace(0, len(data1) - 1, len(data1))

fig, ax = plt.subplots()
ax.plot(x, data1, label="Brute-force Metropolis")
ax.plot(x, data2, label="Metropolis-Hastings")
ax.plot(x, data3, label="Gibbs' sampling")

ax.axhline(asymptote, linestyle='--', color='r', label="Exact")
plt.xlabel("Iteration")
plt.ylabel("Energy")
#plt.axis([-10, 210, 1.5, 6])
#plt.savefig('figure.png')

plt.xlabel(u'Iteration', fontname = 'Times New Roman', size = 14)
plt.ylabel(u'Energy [a.u.]', fontname = 'Times New Roman', size = 14)


#legend = plt.legend(loc='upper right', shadow=False, fontsize='large')
axins = zoomed_inset_axes(ax, 10, loc=10)

axins.axhline(asymptote, linestyle='--', color='r', label="Exact")
axins.plot(x, data1)
axins.plot(x, data2)
axins.plot(x, data3)

x1, x2, y1, y2 = 480, 501, 0.495, 0.505

axins.set_xlim(x1, x2)
axins.set_ylim(y1, y2)

mark_inset(ax, axins, loc1=1, loc2=3, fc="none", ec="0.5")
#axins.grid()
ax.grid()
ax.legend(loc='upper right')

plt.show()
