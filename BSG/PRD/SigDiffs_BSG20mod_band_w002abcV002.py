import matplotlib.pylab as plt
import matplotlib.lines as mlines
import numpy as np
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes 
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

#proton-proton:

x1,y1,ery1 = np.loadtxt("sigma_SD_pp_TOTEM_Cartiglia2013.dat", usecols=(0,1,2), unpack=True)

x2,y2,ery2 = np.loadtxt("sigma_SD_pp_CMS_Khachatryan2015.dat", usecols=(0,1,2), unpack=True)

x3,y3,ery3 = np.loadtxt("sigma_SD_pp_ALICE_Abelev2013.dat", usecols=(0,1,2), unpack=True)

#antiproton-proton:

x4,y4,ery4 = np.loadtxt("sigma_SD_ppbar_UA5_Ansorge1986_Alner1987.dat", usecols=(0,1,2), unpack=True)

x5,y5,ery5 = np.loadtxt("sigma_SD_ppbar_CDF_Abe1994.dat", usecols=(0,1,2), unpack=True)

x6,y6,ery6 = np.loadtxt("sigma_SD_ppbar_E710_Amos1990_Amos1993.dat", usecols=(0,1,2), unpack=True)

x7,y7,ery7 = np.loadtxt("sigma_SD_ppbar_UA4_Bernard1987.dat", usecols=(0,1,2), unpack=True)

x8,y8,ery8 = np.loadtxt("sigma_SD_pp_ISR_Armitage1982.dat", usecols=(0,1,2), unpack=True)

#Sigs Diffs

xf11,yf11 = np.loadtxt("SigSDpp_BSG20mod_Diffs_w002a.dat", usecols=(0,1), unpack=True)

xf12,yf12 = np.loadtxt("SigSDpp_BSG20mod_Diffs_w002b.dat", usecols=(0,1), unpack=True)

xf13,yf13 = np.loadtxt("SigSDpp_BSG20mod_Diffs_w002c.dat", usecols=(0,1), unpack=True)

xf21,yf21 = np.loadtxt("SigSDpbp_BSG20mod_Diffs_w002a.dat", usecols=(0,1), unpack=True)

xf22,yf22 = np.loadtxt("SigSDpbp_BSG20mod_Diffs_w002b.dat", usecols=(0,1), unpack=True)

xf23,yf23 = np.loadtxt("SigSDpbp_BSG20mod_Diffs_w002c.dat", usecols=(0,1), unpack=True)

xf31,yf31 = np.loadtxt("SigDDpp_BSG20mod_Diffs_w002a.dat", usecols=(0,1), unpack=True)

xf32,yf32 = np.loadtxt("SigDDpp_BSG20mod_Diffs_w002b.dat", usecols=(0,1), unpack=True)

xf33,yf33 = np.loadtxt("SigDDpp_BSG20mod_Diffs_w002c.dat", usecols=(0,1), unpack=True)

xf41,yf41 = np.loadtxt("SigDDpbp_BSG20mod_Diffs_w002a.dat", usecols=(0,1), unpack=True)

xf42,yf42 = np.loadtxt("SigDDpbp_BSG20mod_Diffs_w002b.dat", usecols=(0,1), unpack=True)

xf43,yf43 = np.loadtxt("SigDDpbp_BSG20mod_Diffs_w002c.dat", usecols=(0,1), unpack=True)


font = {'family': 'serif','color':  'black', 'weight': 'normal', 'size': 16,}

plt.subplots_adjust(left=0.1, bottom=0.1, right=0.98, top=0.95, wspace=0, hspace=0)

plt.semilogx()
plt.grid(False)
plt.xlim(5,3e4)
plt.ylim(0,25)

p_totem = plt.errorbar(x1,y1,yerr=ery1,fmt='r*',ecolor='r', markeredgecolor='red', capthick=0.5, label=r"TOTEM for $3.4$ GeV $< M_{SD} < 1100$ GeV")
p_cms = plt.errorbar(x2,y2,yerr=ery2,fmt='rs',ecolor='r', markeredgecolor='red', capthick=0.5, label=r"CMS for $M_{X}^{2}/s < 0.05$ or $M_{Y}^{2}/s < 0.05$")
p_alice = plt.errorbar(x3,y3,yerr=ery3,fmt='r^',ecolor='r', markeredgecolor='red', capthick=0.5, label=r"ALICE for $M_{X} < 200$ GeV")
p_isr = plt.errorbar(x8,y8,yerr=ery8,fmt='ko',ecolor='k', capthick=0.5, label=r"ISR for $M_{X}^{2}/s < 0.05$")
p_ua5 = plt.errorbar(x4,y4,yerr=ery4,fmt='w^',ecolor='b', markeredgecolor='blue', capthick=0.5, label=r"UA5 for $M_{X}^{2}/s < 0.05$")
p_cdf = plt.errorbar(x5,y5,yerr=ery5,fmt='wv',ecolor='b', markeredgecolor='blue', capthick=0.5, label=r"CDF for $M_{X}^{2}/s < 0.2$")
p_e710 = plt.errorbar(x6,y6,yerr=ery6,fmt='wo',ecolor='b', markeredgecolor='blue', capthick=0.5, label=r"E710 for $M_{X}^{2}/s < 0.05$")
p_ua4 = plt.errorbar(x7,y7,yerr=ery7,fmt='ws',ecolor='b', markeredgecolor='blue', capthick=0.5, label=r"UA4 for $M_{X}^{2}/s < 0.05$")


plt.fill_between(xf11,yf11,yf13,color='c')
plt.fill_between(xf21,yf21,yf23,color='c')
plt.fill_between(xf31,yf31,yf33,color='c')
plt.fill_between(xf41,yf41,yf43,color='c')

plt.plot(xf12,yf12, 'r-', linewidth=1.25)
plt.plot(xf22,yf22, 'k--', linewidth=1.25)
plt.plot(xf32,yf32, 'r-', linewidth=1.25)
plt.plot(xf42,yf42, 'k--', linewidth=1.25)

graph_data = [p_totem,p_cms,p_alice,p_isr,p_ua5,p_cdf,p_e710,p_ua4]
graph_data_label = ["TOTEM for 3.4 GeV $< M_{SD} <$ 1100 GeV",
  "CMS for $M_{X}^{2}/s < 0.05$ or $M_{Y}^{2}/s < 0.05$", "ALICE for $M_{X} < 200$ GeV",
  "ISR for $M_{X}^{2}/s < 0.05$",
  "UA5 for $M_{X}^{2}/s < 0.05$", "CDF for $M_{X}^{2}/s < 0.2$",
  "E710 for $M_{X}^{2}/s < 0.05$", "UA4 for $M_{X}^{2}/s < 0.05$"]

handle_curve_sd_pp = mlines.Line2D([],[],color='red',linestyle='solid', linewidth=1.25, label=r'$pp$')
handle_curve_sd_ppbar = mlines.Line2D([],[],color='black',linestyle='dashed', linewidth=1.25, label=r'$\bar{p}p$')

artists=[handle_curve_sd_pp,handle_curve_sd_ppbar]#,handle_curve_dd_pp,handle_curve_dd_ppbar]

leg1 = plt.legend(graph_data,graph_data_label,loc=2, ncol=1, shadow=False, fancybox=False, frameon=False, numpoints=1)
leg2 = plt.legend(handles=artists,loc=1, ncol=1, shadow=False, fancybox=False, frameon=False, numpoints=1) # this removes leg1 from the axes.
plt.gca().add_artist(leg1) # add leg1 as a separate artist to the axes

plt.text(8,7,r'$\sigma_{SD}$',fontdict=font)
plt.text(8,2,r'$\sigma_{DD}$',fontdict=font)

plt.xlabel(r"$\sqrt{s}$ [GeV]", fontdict=font)
plt.ylabel(r"$\sigma_{SD}(s),\,\sigma_{DD}(s)$ [mb]", fontdict=font)

plt.savefig("SigDiffs_BSG20mod_band_w002abcV002.eps")
#plt.show()

