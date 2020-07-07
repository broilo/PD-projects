import matplotlib.pylab as plt
import numpy as np
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes 
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

px1,py1 = np.loadtxt("XSecEikPP_BSG20orig_w002b.dat", usecols=(0,1), unpack=True)

#px2,py2 = np.loadtxt("XSecEikPBP_BSG20orig_w002b.dat", usecols=(0,1), unpack=True)

px3,py3 = np.loadtxt("XSecEikPP_BSG20orig_w002b_CTEQ6L.dat", usecols=(0,1), unpack=True)

#px4,py4 = np.loadtxt("XSecEikPBP_BSG20orig_w002b_CTEQ6L.dat", usecols=(0,1), unpack=True)

px5,py5 = np.loadtxt("XSecEikPP_BSG20orig_w002b_MMHT.dat", usecols=(0,1), unpack=True)

#px6,py6 = np.loadtxt("XSecEikPBP_BSG20orig_w002b_MMHT.dat", usecols=(0,1), unpack=True)

font = {'family': 'serif','color':  'black', 'weight': 'normal', 'size': 16,}

plt.subplots_adjust(left=0.1, bottom=0.1, right=0.98, top=0.95, wspace=0, hspace=0)

plt.semilogx()
plt.grid(False)
plt.xlim(5,30000)
plt.ylim(5e1,8e2)

plt.plot(px1,py1, 'k-', linewidth=1.25,label=r"CT14",zorder=3)
#plt.plot(px2,py2, 'k-', linewidth=1.25,zorder=2)

plt.plot(px3,py3, 'r--', linewidth=1.25,label=r"CTEQ6L",zorder=3)
#plt.plot(px4,py4, 'r--', linewidth=1.25,zorder=2)

plt.plot(px5,py5, 'b:', linewidth=1.25,label=r"MMHT",zorder=3)
#plt.plot(px6,py6, 'b:', linewidth=1.25,zorder=2)

plt.xlabel(r"$\sqrt{s}$ [GeV]", fontdict=font)
plt.ylabel(r"$\sigma_{eik}(s)$ [mb]", fontdict=font)

leg = plt.legend(loc=4, ncol=1, shadow=False, fancybox=False, frameon=False, numpoints=1)
leg.get_frame().set_alpha(0.5)

plt.savefig("XSecEik_BSG20orig_w002bPDFs.eps")
#plt.show()

