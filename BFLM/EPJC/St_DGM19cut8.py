import matplotlib.pylab as plt
import numpy as np
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes 
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

data_0,data_1,data_2 = np.loadtxt("paw_StppTAc10GeV.dat", usecols=(0,1,2), unpack=True)

data_3,data_4,data_5 = np.loadtxt("paw_Stpac10GeV.dat", usecols=(0,1,2), unpack=True)

data_6,data_7 = np.loadtxt("Stpp_minDGM18cCTEQ6Llo_sStRh13.dat", usecols=(0,1), unpack=True)
data_8,data_9 = np.loadtxt("Stpa_minDGM18cCTEQ6Llo_sStRh13.dat", usecols=(0,1), unpack=True)

data_10,data_11 = np.loadtxt("Stpp_minDGM18cCT14lo_sStRh13.dat", usecols=(0,1), unpack=True)
data_12,data_13 = np.loadtxt("Stpa_minDGM18cCT14lo_sStRh13.dat", usecols=(0,1), unpack=True)

data_14,data_15 = np.loadtxt("Stpp_minDGM18cMMHTlo_sStRh13.dat", usecols=(0,1), unpack=True)
data_16,data_17 = np.loadtxt("Stpa_minDGM18cMMHTlo_sStRh13.dat", usecols=(0,1), unpack=True)


error1=data_2
error2=data_5

x1=data_6
y1=data_7
x2=data_8
y2=data_9

x3=data_10
y3=data_11
x4=data_12
y4=data_13

x5=data_14
y5=data_15
x6=data_16
y6=data_17

font = {'family': 'serif','color':  'black', 'weight': 'normal', 'size': 16,}
#'color':  'darkred', 'weight': 'normal', 'size': 16,}

plt.subplots_adjust(left=0.1, bottom=0.12, right=0.98, top=0.95, wspace=0, hspace=0)

plt.semilogx()
#plt.grid(True)
plt.xlim(9,20000)
plt.ylim(30,120)
plt.title(r"$\sqrt{s}=8$ TeV", fontdict=font)
#plt.text(2000, 60, 'MSTW', fontdict=font)
#plt.text(2000, 55, '$Q_{min}=1.3$ GeV', fontdict=font)
#plt.text(3000, 52, '$Q_{min}=1.3$ GeV', fontdict=font)
#plt.scatter(data_0,data_1,marker='o',s=10,color='red')
plt.errorbar(data_0,data_1,yerr=error1,fmt='ko',ecolor='k', capthick=0.5)
#plt.scatter(data_2,data_3)
plt.errorbar(data_3,data_4,yerr=error2,fmt='wo',ecolor='k', capthick=0.5)

plt.plot(data_6,data_7, 'k-', linewidth=1.5, label='CTEQ6L')
plt.plot(data_8,data_9, 'k-', linewidth=1.5)

plt.plot(data_10,data_11, 'b--', linewidth=1.5, label='CTEQ6L1')
plt.plot(data_12,data_13, 'b--', linewidth=1.5)

plt.plot(data_14,data_15, 'r:', linewidth=1.7, label='MMHT')
plt.plot(data_16,data_17, 'r:', linewidth=1.7)


plt.xlabel(r"$\sqrt{s}$ [GeV]", fontsize=20)
plt.ylabel(r"$\sigma_{tot}(s)$ [mb]", fontsize=20)

leg = plt.legend(loc=4, ncol=1, shadow=True, fancybox=True, frameon=True, numpoints=1)
leg.get_frame().set_alpha(0.5)

###########################################
sub_axes = plt.axes([.2, .55, .35, .35])
#plt.grid(True)
plt.semilogx()
plt.xlim(5000,15000)
plt.ylim(90,120)
plt.errorbar(data_0,data_1,yerr=error1,fmt='ko',ecolor='k', capthick=0.5)
plt.errorbar(data_3,data_4,yerr=error2,fmt='wo',ecolor='k', capthick=0.5)
plt.plot(data_6,data_7, 'k-', linewidth=1.5)
plt.plot(data_8,data_9, 'k-', linewidth=1.5)
plt.plot(data_10,data_11, 'b--', linewidth=1.5)
plt.plot(data_12,data_13, 'b--', linewidth=1.5)
plt.plot(data_14,data_15, 'r:', linewidth=1.5)
plt.plot(data_16,data_17, 'r:', linewidth=1.5)
###########################################

plt.savefig("St_DGM19cut8.eps")
plt.show()
