import matplotlib.pylab as plt
import numpy as np
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes 
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

px1,py1 = np.loadtxt("XSecSoftPP_BSG20orig_w002b.dat", usecols=(0,1), unpack=True)

#px2,py2 = np.loadtxt("XSecSoftPBP_BSG20orig_w002b.dat", usecols=(0,1), unpack=True)

px3,py3 = np.loadtxt("XSecSoftPP_BSG20orig_w002b_NNPDF31.dat", usecols=(0,1), unpack=True)

#px4,py4 = np.loadtxt("XSecSoftPBP_BSG20orig_w002b_NNPDF31.dat", usecols=(0,1), unpack=True)

px5,py5 = np.loadtxt("XSecSoftPP_BSG20orig_w002b_MMHT.dat", usecols=(0,1), unpack=True)

#px6,py6 = np.loadtxt("XSecSoftPBP_BSG20orig_w002b_MMHT.dat", usecols=(0,1), unpack=True)

px7,py7 = np.loadtxt("XSecHard_BSG20orig_w002b.dat", usecols=(0,1), unpack=True)

px8,py8 = np.loadtxt("XSecHard_BSG20orig_w002b_NNPDF31.dat", usecols=(0,1), unpack=True)

px9,py9 = np.loadtxt("XSecHard_BSG20orig_w002b_MMHT.dat", usecols=(0,1), unpack=True)

font = {'family': 'serif','color':  'black', 'weight': 'normal', 'size': 16,}

plt.subplots_adjust(left=0.1, bottom=0.1, right=0.98, top=0.95, wspace=0, hspace=0)

plt.semilogx()
plt.grid(False)
plt.xlim(5,30000)
plt.ylim(0,2.0e2)

plt.plot(px1,py1, 'k-', linewidth=1.25,label="CT14",zorder=3)
#plt.plot(px2,py2, 'k-', linewidth=1.25,zorder=2)

plt.plot(px3,py3, 'r--', linewidth=1.25,label="NNPDF31",zorder=3)
#plt.plot(px4,py4, 'r--', linewidth=1.25,zorder=2)

plt.plot(px5,py5, 'b:', linewidth=1.25,label="MMHT",zorder=3)
#plt.plot(px6,py6, 'b:', linewidth=1.25,zorder=2)

plt.plot(px7,py7, 'k-', linewidth=1.25)
plt.plot(px8,py8, 'r--', linewidth=1.25)
plt.plot(px9,py9, 'b:', linewidth=1.25)

plt.xlabel(r"$\sqrt{s}$ [GeV]", fontdict=font)
plt.ylabel(r"$\sigma(s)$ [mb]", fontdict=font)

leg = plt.legend(loc=4, ncol=1, shadow=False, fancybox=False, frameon=False, numpoints=1)
leg.get_frame().set_alpha(0.5)

plt.savefig("XSecSoftHard_BSG20orig_w002bPDFsV2_v3.eps")
#plt.show()

