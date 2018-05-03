import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes, mark_inset

asymptote = 0.5

data1 = np.loadtxt("../data/BF_1P_N_1_MC_5000000_eta_005.txt")
data2 = np.loadtxt("../data/BF_1P_N_1_MC_5000000_eta_001.txt")
data3 = np.loadtxt("../data/BF_1P_N_1_MC_5000000_eta_05.txt")
#data5 = np.loadtxt("../data/BF_1P_N_1_MC_5000000_eta_1.txt")
#data6 = np.loadtxt("../data/BF_1P_N_1_MC_5000000_eta_16.txt")

x = np.linspace(0, len(data1) - 1, len(data1))

fig, ax = plt.subplots()
ax.axhline(asymptote, linestyle='--', color='r', label="Exact")
ax.plot(x, data1, label = '$\eta=0.05$')
ax.plot(x, data2, label = '$\eta=0.01$')
ax.plot(x, data3, label = '$\eta=0.5$')
#ax.plot(x, data5, label = '$\eta=1.0$')
#ax.plot(x, data6, label = '$\eta=1.5$')

plt.xlabel(u'Iteration', fontname = 'Times New Roman', size = 14)
plt.ylabel(u'Energy [a.u.]', fontname = 'Times New Roman', size = 14)


#legend = plt.legend(loc='upper right', shadow=False, fontsize='large')
axins = zoomed_inset_axes(ax, 600, loc=9)

axins.axhline(asymptote, linestyle='--', color='r', label="Exact")
axins.plot(x, data1, label = '$\eta=0.05$')
axins.plot(x, data2, label = '$\eta=0.01$')
axins.plot(x, data3, label = '$\eta=0.5$')
#axins.plot(x, data5, label = '$\eta=1.0$')
#axins.plot(x, data6, label = '$\eta=1.5$')

x1, x2, y1, y2 = 98.96, 99.01, 0.49999, 0.50001

axins.set_xlim(x1, x2)
axins.set_ylim(y1, y2)

mark_inset(ax, axins, loc1=1, loc2=3, fc="none", ec="0.5")
#axins.grid()
ax.grid()
ax.legend(loc='upper right')

plt.show()
