import matplotlib.pylab as plt
import numpy as np
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes 
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

########################################################################
#### Conjunto de dados
#
data_1,data_2,data_3 = np.loadtxt("sig_sumij_real_CTEQ6Llo_m400q13V001.dat", usecols=(0,1,2), unpack=True)

x1=data_1
y1=data_2
error1=data_3

data_4,data_5,data_6 = np.loadtxt("sig_sumij_real_CT14lo_m400q13V001.dat", usecols=(0,1,2), unpack=True)

x2=data_4
y2=data_5
error2=data_6

data_7,data_8,data_9 = np.loadtxt("sig_sumij_real_MMHTlo_m400q13V001.dat", usecols=(0,1,2), unpack=True)

x3=data_7
y3=data_8
error3=data_9
#
########################################################################
########################################################################
#### Fits
#
data_11,data_12 = np.loadtxt("gerReSigQCDder_CTEQ6Llo.dat", usecols=(0,1), unpack=True)
x11=data_11
y11=data_12

data_13,data_14 = np.loadtxt("gerImSigQCDder_CTEQ6Llo.dat", usecols=(0,1), unpack=True)
x12=data_13
y12=data_14

data_41,data_42 = np.loadtxt("gerReSigQCDder_CT14lo.dat", usecols=(0,1), unpack=True)
x21=data_41
y21=data_42

data_43,data_44 = np.loadtxt("gerImSigQCDder_CT14lo.dat", usecols=(0,1), unpack=True)
x22=data_43
y22=data_44

data_71,data_72 = np.loadtxt("gerReSigQCDder_MMHTlo.dat", usecols=(0,1), unpack=True)
x71=data_71
y71=data_72

data_73,data_74 = np.loadtxt("gerImSigQCDder_MMHTlo.dat", usecols=(0,1), unpack=True)
x72=data_73
y72=data_74
#
########################################################################

font = {'family': 'serif','color':  'black', 'weight': 'normal', 'size': 16,}
#'color':  'darkred', 'weight': 'normal', 'size': 16,}

plt.subplots_adjust(left=0.13, bottom=0.12, right=0.95, top=0.95, wspace=0, hspace=0)

plt.semilogx()
plt.grid(True)
plt.xlim(8,1000000)
plt.ylim(0,15000)
plt.title('Leading order', fontdict=font)
plt.text(20, 7000, 'Real parametrization', fontdict=font)
plt.text(20, 6000, '$Q_{min}=1.3$ GeV', fontdict=font)
plt.text(20, 5000, '$m_{g}=400$ MeV', fontdict=font)
########################################################################
#### Conjunto de dados
#
plt.errorbar(x1,y1,yerr=error1,fmt='k.',ecolor='k', capthick=1.8)
plt.errorbar(x2,y2,yerr=error2,fmt='b.',ecolor='b', capthick=1.8)
plt.errorbar(x3,y3,yerr=error3,fmt='r.',ecolor='r', capthick=1.8)
#
########################################################################
#### Fits
#
plt.plot(x11,y11, 'k-', linewidth=1.5, label=r'(CTEQ6L) $Re\,\sigma_{QCD}(s)$ Parametrized')
plt.plot(x12,y12, 'k--', linewidth=1.5, label=r'(CTEQ6L) $-Im\,\sigma_{QCD}(s)$ RDD')

plt.plot(x21,y21, 'b-', linewidth=1.5, label=r'(CT14) $Re\,\sigma_{QCD}(s)$ Parametrized')
plt.plot(x22,y22, 'b--', linewidth=1.5, label=r'(CT14) $-Im\,\sigma_{QCD}(s)$ RDD')

plt.plot(x71,y71, 'r-', linewidth=1.5, label=r'(MMHT) $Re\,\sigma_{QCD}(s)$ Parametrized')
plt.plot(x72,y72, 'r--', linewidth=1.5, label=r'(MMHT) $-Im\,\sigma_{QCD}(s)$ RDD')
#
########################################################################

plt.xlabel(r"$\sqrt{s} \,\,\, [GeV]$", fontsize=20)
plt.ylabel(r"$\sigma_{QCD}(s)\, [GeV^{-2}]$", fontsize=20)


leg = plt.legend(loc=2, ncol=1, shadow=True, fancybox=True)
leg.get_frame().set_alpha(0.5)

plt.savefig("SigsQCDder.eps")
plt.show()
