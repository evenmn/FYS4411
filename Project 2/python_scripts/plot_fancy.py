import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes, mark_inset

asymptote = 0.5

data1 = np.loadtxt("../data/BF_1P_N_1_MC_1000000_eta_001.txt")
data2 = np.loadtxt("../data/BF_1P_N_2_MC_1000000_eta_001.txt")
data3 = np.loadtxt("../data/BF_1P_N_3_MC_1000000_eta_001.txt")
data5 = np.loadtxt("../data/BF_1P_N_5_MC_1000000_eta_001.txt")
data6 = np.loadtxt("../data/BF_1P_N_10_MC_1000000_eta_001.txt")

x = np.linspace(0, len(data1) - 1, len(data1))

fig, ax = plt.subplots()
ax.axhline(asymptote, linestyle='--', color='r', label="Exact")
ax.plot(x, data1, label = 'N=1')
ax.plot(x, data2, label = 'N=2')
ax.plot(x, data3, label = 'N=3')
ax.plot(x, data5, label = 'N=5')
ax.plot(x, data6, label = 'N=10')

plt.xlabel(u'Iteration', fontname = 'Times New Roman', size = 14)
plt.ylabel(u'Energy [a.u.]', fontname = 'Times New Roman', size = 14)

#legend = plt.legend(loc='upper right', shadow=False, fontsize='large')
axins = zoomed_inset_axes(ax, 8, loc=10)

axins.axhline(asymptote, linestyle='--', color='r', label="Exact")
axins.plot(x, data1, label = 'N=1')
axins.plot(x, data2, label = 'N=2')
axins.plot(x, data3, label = 'N=3')
axins.plot(x, data5, label = 'N=5')
axins.plot(x, data6, label = 'N=10')

x1, x2, y1, y2 = 95, 100, 0.47, 0.58

axins.set_xlim(x1, x2)
axins.set_ylim(y1, y2)

mark_inset(ax, axins, loc1=1, loc2=3, fc="none", ec="0.5")

ax.grid()
#axins.grid()
ax.legend(loc='upper right')
#ax.axis([0, 110, 0.8, 0.4])
#plt.savefig('figure.png')

plt.show()
