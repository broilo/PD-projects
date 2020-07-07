import matplotlib.pylab as plt
import numpy as np
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes 
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

x1,y1,erry1 = np.loadtxt("paw_Stpa.dat", usecols=(0,1,2), unpack=True)

x2,y2,erry2 = np.loadtxt("paw_Stpp.dat", usecols=(0,1,2), unpack=True)


px1,py1 = np.loadtxt("Stpp_BSG20mod_w002b.dat", usecols=(0,1), unpack=True)

px2,py2 = np.loadtxt("Stpbp_BSG20mod_w002b.dat", usecols=(0,1), unpack=True)


px3,py3 = np.loadtxt("Stpp_BSG20mod_w002a.dat", usecols=(0,1), unpack=True)

px4,py4 = np.loadtxt("Stpbp_BSG20mod_w002a.dat", usecols=(0,1), unpack=True)


px5,py5 = np.loadtxt("Stpp_BSG20mod_w002a1_5.dat", usecols=(0,1), unpack=True)

px6,py6 = np.loadtxt("Stpbp_BSG20mod_w002a1_5.dat", usecols=(0,1), unpack=True)


px7,py7 = np.loadtxt("Stpp_BSG20mod_w002a1_6.dat", usecols=(0,1), unpack=True)

px8,py8 = np.loadtxt("Stpbp_BSG20mod_w002a1_6.dat", usecols=(0,1), unpack=True)


px9,py9 = np.loadtxt("Stpp_BSG20mod_w002c.dat", usecols=(0,1), unpack=True)

px10,py10 = np.loadtxt("Stpbp_BSG20mod_w002c.dat", usecols=(0,1), unpack=True)


px11,py11 = np.loadtxt("Stpp_BSG20mod_w002c1_8.dat", usecols=(0,1), unpack=True)

px12,py12 = np.loadtxt("Stpbp_BSG20mod_w002c1_8.dat", usecols=(0,1), unpack=True)


px13,py13 = np.loadtxt("Stpp_BSG20mod_w002c1_9.dat", usecols=(0,1), unpack=True)

px14,py14 = np.loadtxt("Stpbp_BSG20mod_w002c1_9.dat", usecols=(0,1), unpack=True)


font = {'family': 'serif','color':  'black', 'weight': 'normal', 'size': 16,}

plt.subplots_adjust(left=0.1, bottom=0.1, right=0.98, top=0.95, wspace=0, hspace=0)

plt.semilogx()
plt.grid(False)
plt.xlim(5,30000)
plt.ylim(30,130)


plt.errorbar(x1,y1,yerr=erry1,fmt='wo',ecolor='k', capthick=0.5, label=r"$\bar{p}p$")
plt.errorbar(x2,y2,yerr=erry2,fmt='ko',ecolor='k', capthick=0.5, label=r"$pp$")


plt.plot(px3,py3, 'b-', linewidth=0.5,label=r"$w=1.4$")
plt.plot(px4,py4, 'b--', linewidth=0.5)

plt.plot(px5,py5, 'g-', linewidth=0.5,label=r"$w=1.5$")
plt.plot(px6,py6, 'g--', linewidth=0.5)

plt.plot(px7,py7, 'r-', linewidth=0.5,label=r"$w=1.6$")
plt.plot(px8,py8, 'r--', linewidth=0.5)

plt.plot(px1,py1, 'k-', linewidth=0.5,label=r"$w=1.72$")
plt.plot(px2,py2, 'k--', linewidth=0.5)

plt.plot(px11,py11, 'c-', linewidth=0.5,label=r"$w=1.8$")
plt.plot(px12,py12, 'c--', linewidth=0.5)

plt.plot(px13,py13, 'm-', linewidth=0.5,label=r"$w=1.9$")
plt.plot(px14,py14, 'm--', linewidth=0.5)

plt.plot(px9,py9, 'y-', linewidth=0.5,label=r"$w=2.0$")
plt.plot(px10,py10, 'y--', linewidth=0.5)

plt.xlabel(r"$\sqrt{s}$ [GeV]", fontdict=font)
plt.ylabel(r"$\sigma_{tot}(s)$ [mb]", fontdict=font)

leg = plt.legend(loc=4, ncol=1, shadow=False, fancybox=False, frameon=False, numpoints=1)
leg.get_frame().set_alpha(0.5)

###########################################
sub_axes = plt.axes([.2, .55, .35, .35])
plt.grid(False)
plt.semilogx()
plt.xlim(5000,15000)
plt.ylim(90,120)

plt.errorbar(x1,y1,yerr=erry1,fmt='wo',ecolor='k', capthick=0.5)
plt.errorbar(x2,y2,yerr=erry2,fmt='ko',ecolor='k', capthick=0.5)

plt.plot(px3,py3, 'b-', linewidth=0.5,label=r"$w=1.4$")
plt.plot(px4,py4, 'b--', linewidth=0.5)

plt.plot(px5,py5, 'g-', linewidth=0.5,label=r"$w=1.5$")
plt.plot(px6,py6, 'g--', linewidth=0.5)

plt.plot(px7,py7, 'r-', linewidth=0.5,label=r"$w=1.6$")
plt.plot(px8,py8, 'r--', linewidth=0.5)

plt.plot(px1,py1, 'k-', linewidth=0.5,label=r"$w=1.72$")
plt.plot(px2,py2, 'k--', linewidth=0.5)

plt.plot(px11,py11, 'c-', linewidth=0.5,label=r"$w=1.8$")
plt.plot(px12,py12, 'c--', linewidth=0.5)

plt.plot(px13,py13, 'm-', linewidth=0.5,label=r"$w=1.9$")
plt.plot(px14,py14, 'm--', linewidth=0.5)

plt.plot(px9,py9, 'y-', linewidth=0.5,label=r"$w=2.0$")
plt.plot(px10,py10, 'y--', linewidth=0.5)

plt.plot(px1,py1, 'k-', linewidth=1.0)
plt.plot(px2,py2, 'k--', linewidth=1.0)

###########################################

plt.savefig("St_BSG20mod_w002s.eps")
#plt.show()

