import matplotlib.pylab as plt
import numpy as np
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes 
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

"""
Essa figura eu usei o programa BSG20orig_w002b.C usando os ajustes de BGS20mod do paper
"""

px1,py1 = np.loadtxt("XSecEikPP_BSG20orig_w002b.dat", usecols=(0,1), unpack=True)

px2,py2 = np.loadtxt("XSecEikPBP_BSG20orig_w002b.dat", usecols=(0,1), unpack=True)

px3,py3 = np.loadtxt("XSecSoftPP_BSG20orig_w002b.dat", usecols=(0,1), unpack=True)

px4,py4 = np.loadtxt("XSecSoftPBP_BSG20orig_w002b.dat", usecols=(0,1), unpack=True)

px5,py5 = np.loadtxt("XSecHard_BSG20orig_w002b.dat", usecols=(0,1), unpack=True)

font = {'family': 'serif','color':  'black', 'weight': 'normal', 'size': 16,}

plt.subplots_adjust(left=0.1, bottom=0.1, right=0.98, top=0.95, wspace=0, hspace=0)

plt.semilogx()
plt.grid(False)
plt.xlim(5,30000)
plt.ylim(0,800)


plt.plot(px1,py1, 'r-', linewidth=1.25,label=r"$\sigma^{pp}_{eik}$",zorder=2)
plt.plot(px2,py2, 'k--', linewidth=1.25,label=r"$\sigma^{\bar{p}p}_{eik}$",zorder=3)

plt.plot(px3,py3, 'r:', linewidth=1.25,label=r"$\sigma^{pp}_{soft}$",zorder=3)
plt.plot(px4,py4, 'k:', linewidth=1.25,label=r"$\sigma^{\bar{p}p}_{soft}$",zorder=3)

plt.plot(px5,py5, 'b-', linewidth=1.25,label=r"$\sigma_{pQCD}$",zorder=3)

plt.xlabel(r"$\sqrt{s}$ [GeV]", fontdict=font)
plt.ylabel(r"$\sigma(s)$ [mb]", fontdict=font)

leg = plt.legend(loc=2, ncol=1, shadow=False, fancybox=False, frameon=False, numpoints=1)
leg.get_frame().set_alpha(0.5)



plt.savefig("XSec_eik_BSG20mod_changes.eps")
#plt.show()

