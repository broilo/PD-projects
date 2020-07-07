import matplotlib.pylab as plt
import numpy as np
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes 
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

data_0,data_1,data_2 = np.loadtxt("paw_Rhpp.dat", usecols=(0,1,2), unpack=True)

data_3,data_4,data_5 = np.loadtxt("paw_Rhpa.dat", usecols=(0,1,2), unpack=True)

data_6,data_7 = np.loadtxt("RhPP_BSG20orig_w002b_F77.dat", usecols=(0,1), unpack=True)
data_8,data_9 = np.loadtxt("RhPBP_BSG20orig_w002b_F77.dat", usecols=(0,1), unpack=True)

error1=data_2
error2=data_5

x1=data_6
y1=data_7
x2=data_8
y2=data_9

font = {'family': 'serif','color':  'black', 'weight': 'normal', 'size': 16,}
#'color':  'darkred', 'weight': 'normal', 'size': 16,}

plt.subplots_adjust(left=0.11, bottom=0.12, right=0.98, top=0.95, wspace=0, hspace=0)

plt.semilogx()
#plt.grid(True)
plt.xlim(5,30000)
plt.ylim(-0.3,0.3)
plt.title('Born-level', fontdict=font)
plt.errorbar(data_0,data_1,yerr=error1,fmt='ko',ecolor='k', capthick=0.5, label=r'$pp$')
plt.errorbar(data_3,data_4,yerr=error2,fmt='wo',ecolor='k', capthick=0.5, label=r'$\bar{p}p$')

plt.plot(data_8,data_9, 'r-', linewidth=1.5, label=r'$\bar{p}p$')
plt.plot(data_6,data_7, 'b-', linewidth=1.5, label=r'$pp$')

plt.xlabel(r"$\sqrt{s}$ [GeV]", fontdict=font)
plt.ylabel(r"$\rho(s)$", fontdict=font)

leg = plt.legend(loc=2, ncol=1, shadow=False, fancybox=False, frameon=False, numpoints=1)
leg.get_frame().set_alpha(0.5)

###########################################
sub_axes = plt.axes([.6, .2, .35, .35])
#plt.grid(True)
plt.semilogx()
plt.xlim(5000,15000)
plt.ylim(0.07,0.16)
plt.errorbar(data_0,data_1,yerr=error1,fmt='ko',ecolor='k', capthick=0.5)
plt.errorbar(data_3,data_4,yerr=error2,fmt='wo',ecolor='k', capthick=0.5)
plt.plot(data_8,data_9, 'r-', linewidth=1.5)
plt.plot(data_6,data_7, 'b-', linewidth=1.5)
###########################################

plt.savefig("Rh_BSG20orig_w002b_F77.eps")
#plt.show()
