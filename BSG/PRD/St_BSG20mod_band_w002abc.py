import matplotlib.pylab as plt
import numpy as np
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes 
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

x1,y1,erry1 = np.loadtxt("paw_Stpa.dat", usecols=(0,1,2), unpack=True)

x2,y2,erry2 = np.loadtxt("paw_Stpp_semLHC.dat", usecols=(0,1,2), unpack=True)

x3,y3,erry3 = np.loadtxt("paw_StATLAS.dat", usecols=(0,1,2), unpack=True)

x4,y4,erry4 = np.loadtxt("paw_StTOTEM.dat", usecols=(0,1,2), unpack=True)

px1,py1 = np.loadtxt("Stpp_BSG20mod_w002b.dat", usecols=(0,1), unpack=True)

px2,py2 = np.loadtxt("Stpbp_BSG20mod_w002b.dat", usecols=(0,1), unpack=True)

px3,py3 = np.loadtxt("Stpp_BSG20mod_w002a.dat", usecols=(0,1), unpack=True)

px4,py4 = np.loadtxt("Stpp_BSG20mod_w002c.dat", usecols=(0,1), unpack=True)

px5,py5 = np.loadtxt("Stpbp_BSG20mod_w002a.dat", usecols=(0,1), unpack=True)

px6,py6 = np.loadtxt("Stpbp_BSG20mod_w002c.dat", usecols=(0,1), unpack=True)

font = {'family': 'serif','color':  'black', 'weight': 'normal', 'size': 16,}

plt.subplots_adjust(left=0.1, bottom=0.1, right=0.98, top=0.95, wspace=0, hspace=0)

plt.semilogx()
plt.grid(False)
plt.xlim(5,30000)
plt.ylim(30,130)


plt.errorbar(x1,y1,yerr=erry1,fmt='wo',ecolor='k', capthick=0.5, label=r"$\bar{p}p$",zorder=1)
plt.errorbar(x2,y2,yerr=erry2,fmt='ko',ecolor='k', capthick=0.5, label=r"$pp$ below LHC",zorder=1)
plt.errorbar(x3,y3,yerr=erry3,fmt='bs',ecolor='b', markeredgecolor='blue', capthick=0.5, label="ATLAS",zorder=1)
plt.errorbar(x4,y4,yerr=erry4,fmt='r^',ecolor='r', markeredgecolor='red', capthick=0.5, label="TOTEM",zorder=1)


plt.fill_between(px3,py3,py4,color='c',zorder=0)
plt.fill_between(px5,py5,py6,color='c',zorder=0)

plt.plot(px1,py1, 'r-', linewidth=1.25,label=r"$pp$",zorder=3)
plt.plot(px2,py2, 'k--', linewidth=1.25,label=r"$\bar{p}p$",zorder=2)

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
plt.errorbar(x3,y3,yerr=erry3,fmt='bs',ecolor='b', markeredgecolor='blue', capthick=0.5)
plt.errorbar(x4,y4,yerr=erry4,fmt='r^',ecolor='r', markeredgecolor='red', capthick=0.5)

plt.fill_between(px3,py3,py4,color='c')
plt.fill_between(px5,py5,py6,color='c')

plt.plot(px1,py1, 'r-', linewidth=1.0)
plt.plot(px2,py2, 'k--', linewidth=1.0)

###########################################

plt.savefig("St_BSG20mod_band_w002abc.eps")
#plt.show()

