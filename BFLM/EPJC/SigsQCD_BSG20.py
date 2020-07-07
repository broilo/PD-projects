import matplotlib.pylab as plt
import numpy as np
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes 
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

########################################################################
#### Conjunto de dados
#
x1,y1,ery1 = np.loadtxt("sig_sumij_real_CTEQ6Llo_m400q13V001.dat", usecols=(0,1,2), unpack=True)

x2,y2,ery2 = np.loadtxt("sig_sumij_real_CT14lo_m400q13V001.dat", usecols=(0,1,2), unpack=True)

x3,y3,ery3 = np.loadtxt("sig_sumij_real_MMHTlo_m400q13V001.dat", usecols=(0,1,2), unpack=True)
#
########################################################################
########################################################################
#### Fits
#
xf1,yf1 = np.loadtxt("gerReSigQCDpresc_CTEQ6Llo.dat", usecols=(0,1), unpack=True)

xf2,yf2 = np.loadtxt("gerImSigQCDpresc_CTEQ6Llo.dat", usecols=(0,1), unpack=True)

xf3,yf3 = np.loadtxt("gerReSigQCDpresc_CT14lo.dat", usecols=(0,1), unpack=True)

xf4,yf4 = np.loadtxt("gerImSigQCDpresc_CT14lo.dat", usecols=(0,1), unpack=True)

xf5,yf5 = np.loadtxt("gerReSigQCDpresc_MMHTlo.dat", usecols=(0,1), unpack=True)

xf6,yf6 = np.loadtxt("gerImSigQCDpresc_MMHTlo.dat", usecols=(0,1), unpack=True)
#
########################################################################
font = {'family': 'serif','color':  'black', 'weight': 'normal', 'size': 16,}

plt.subplots_adjust(left=0.13, bottom=0.12, right=0.95, top=0.95, wspace=0, hspace=0)

plt.semilogx()
#plt.grid(True)
plt.xlim(8,100000)
plt.ylim(0,8000)
#plt.title('Leading order', fontdict=font)
#plt.text(20, 5800, 'Complex parametrization', fontdict=font)
plt.text(20, 6500, '$Q_{min}=1.3$ GeV', fontdict=font)
plt.text(20, 6000, '$m_{g}=400$ MeV', fontdict=font)
########################################################################
#### Conjunto de dados
#
#plt.errorbar(x1,y1,yerr=ery1,fmt='k.',ecolor='k', capthick=1.8)
plt.errorbar(x2,y2,yerr=ery2,fmt='r.',ecolor='r', capthick=1.8, label="CT14 simulated points")
#plt.errorbar(x3,y3,yerr=ery3,fmt='r.',ecolor='r', capthick=1.8)
#
########################################################################
#### Fits
#
#plt.plot(xf1,yf1, 'k-', linewidth=1.5, label=r'(CTEQ6L) $Re\,\sigma_{QCD}(s)$ Parametrized')
#plt.plot(xf2,yf2, 'k--', linewidth=1.5, label=r'(CTEQ6L) $-Im\,\sigma_{QCD}(s)$ $s\to-is$')

plt.plot(xf3,yf3, 'k-', linewidth=1.0, label=r'$\sigma_{pQCD}(s)$ Parametrized')
#plt.plot(xf4,yf4, 'k--', linewidth=1.5, label=r'$-Im\,\sigma_{PQCD}(s)$ $s\to-is$')

#plt.plot(xf5,yf5, 'r-', linewidth=1.5, label=r'(MMHT) $Re\,\sigma_{QCD}(s)$ Parametrized')
#plt.plot(xf6,yf6, 'r--', linewidth=1.5, label=r'(MMHT) $-Im\,\sigma_{QCD}(s)$ $s\to-is$')
#
########################################################################

plt.xlabel(r"$\sqrt{s} \,\,\, [GeV]$", fontdict=font)
plt.ylabel(r"$\sigma_{pQCD}(s)\, [GeV^{-2}]$", fontdict=font)

leg = plt.legend(loc=2, ncol=1, shadow=False, fancybox=False, frameon=False, numpoints=1)
leg.get_frame().set_alpha(0.5)

plt.savefig("SigsQCD_BSG20.eps")
#plt.show()
