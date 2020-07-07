import matplotlib.pylab as plt
import numpy as np
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes 
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

data_0,data_1 = np.loadtxt("xgVSx_CTEQ6L_Q10.dat", usecols=(0,1), unpack=True)
data_2,data_3 = np.loadtxt("xgVSx_CT14_Q10.dat", usecols=(0,1), unpack=True)
data_4,data_5 = np.loadtxt("xgVSx_MMHT_Q10.dat", usecols=(0,1), unpack=True)

x1=data_0
y1=0.0025*data_1
x2=data_2
y2=0.0025*data_3
x3=data_4
y3=0.0025*data_5

font = {'family': 'serif','color':  'black', 'weight': 'normal', 'size': 16,}

plt.subplots_adjust(left=0.13, bottom=0.11, right=0.98, top=0.95, wspace=0, hspace=0)

plt.semilogx()
#plt.grid(True)
plt.xlim(0.00001,1)
plt.ylim(0,1.0)

plt.plot(x1,y1, 'k-', linewidth=2.8, label='CTEQ6L')
plt.plot(x2,y2, 'b--', linewidth=2.8, label='CT14')
plt.plot(x3,y3, 'r:', linewidth=2.8, label='MMHT')

plt.title('$Q=10$ GeV', fontdict=font)
plt.xlabel(r"$x$", fontsize=20)
plt.ylabel(r"$xg(x,Q^{2}) [\times 0.0025]$", fontsize=20)

leg = plt.legend(loc=1, ncol=1, shadow=True, fancybox=True, numpoints=1, frameon=True)
leg.get_frame().set_alpha(0.5)

plt.savefig("xgVSx_Q10v2.eps")
plt.show()

