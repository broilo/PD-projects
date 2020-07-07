import matplotlib.pylab as plt
import numpy as np
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes 
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

x1,y1 = np.loadtxt("BSG20mod_chi2dof_w002.dat", usecols=(0,1), unpack=True)

x2,y2 = np.loadtxt("BSG20mod_chi2dof_w002Fit.dat", usecols=(0,1), unpack=True)

font = {'family': 'serif','color':  'black', 'weight': 'normal', 'size': 16,}

plt.subplots_adjust(left=0.1, bottom=0.1, right=0.98, top=0.95, wspace=0, hspace=0)


plt.grid(False)
plt.xlim(1.2,2.2)
plt.ylim(0,170)

plt.errorbar(x1,y1,fmt='ko', capthick=0.5)
plt.errorbar(x2,y2,fmt='rs', capthick=0.5, label=r"(1.72,3.31)")

plt.xlabel(r"$w$", fontdict=font)
plt.ylabel(r"$\chi^{2}/dof$ ", fontdict=font)

leg = plt.legend(loc=4, ncol=1, shadow=False, fancybox=False, frameon=False, numpoints=1)
leg.get_frame().set_alpha(0.5)

###########################################
#sub_axes = plt.axes([.2, .55, .35, .35])
sub_axes = plt.axes([.65, .64, .28, .28])
plt.grid(False)
plt.xlim(1.7,1.75)
plt.ylim(3.2,4)

plt.errorbar(x1,y1,fmt='ko', capthick=0.5)
plt.errorbar(x2,y2,fmt='rs', capthick=0.5)

###########################################

plt.savefig("chi2dof_BSG20mod_w002.eps")
#plt.show()

