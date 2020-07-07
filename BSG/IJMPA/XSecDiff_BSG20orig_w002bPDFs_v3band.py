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

px1,py1 = np.loadtxt("XSecDiffPP_BSG20orig_w002b.dat", usecols=(0,1), unpack=True)

px2,py2 = np.loadtxt("XSecDiffPBP_BSG20orig_w002b.dat", usecols=(0,1), unpack=True)


px3,py3 = np.loadtxt("XSecDiffPP_BSG20orig_w002b_NNPDF31.dat", usecols=(0,1), unpack=True)
px31,py31 = np.loadtxt("XSecDiffPP_BSG20orig_w002b_NNPDF31_member_38.dat", usecols=(0,1), unpack=True)
px32,py32 = np.loadtxt("XSecDiffPP_BSG20orig_w002b_NNPDF31_member_51.dat", usecols=(0,1), unpack=True)

px4,py4 = np.loadtxt("XSecDiffPBP_BSG20orig_w002b_NNPDF31.dat", usecols=(0,1), unpack=True)
px41,py41 = np.loadtxt("XSecDiffPBP_BSG20orig_w002b_NNPDF31_member_38.dat", usecols=(0,1), unpack=True)
px42,py42 = np.loadtxt("XSecDiffPBP_BSG20orig_w002b_NNPDF31_member_51.dat", usecols=(0,1), unpack=True)

px5,py5 = np.loadtxt("XSecDiffPP_BSG20orig_w002b_MMHT.dat", usecols=(0,1), unpack=True)
px51,py51 = np.loadtxt("XSecDiffPP_BSG20orig_w002b_MMHT_member_34.dat", usecols=(0,1), unpack=True)
px52,py52 = np.loadtxt("XSecDiffPP_BSG20orig_w002b_MMHT_member_39.dat", usecols=(0,1), unpack=True)

px6,py6 = np.loadtxt("XSecDiffPBP_BSG20orig_w002b_MMHT.dat", usecols=(0,1), unpack=True)
px61,py61 = np.loadtxt("XSecDiffPBP_BSG20orig_w002b_MMHT_member_34.dat", usecols=(0,1), unpack=True)
px62,py62 = np.loadtxt("XSecDiffPBP_BSG20orig_w002b_MMHT_member_39.dat", usecols=(0,1), unpack=True)


font = {'family': 'serif','color':  'black', 'weight': 'normal', 'size': 16,}

plt.subplots_adjust(left=0.1, bottom=0.1, right=0.98, top=0.95, wspace=0, hspace=0)

plt.semilogx()
plt.grid(False)
plt.xlim(5,3e4)
plt.ylim(0,30)

p_totem = plt.errorbar(x1,y1,yerr=ery1,fmt='r*',ecolor='r', markeredgecolor='red', capthick=0.5, label=r"TOTEM for $3.4$ GeV $< M_{SD} < 1100$ GeV")
p_cms = plt.errorbar(x2,y2,yerr=ery2,fmt='rs',ecolor='r', markeredgecolor='red', capthick=0.5, label=r"CMS for $M_{X}^{2}/s < 0.05$ or $M_{Y}^{2}/s < 0.05$")
p_alice = plt.errorbar(x3,y3,yerr=ery3,fmt='r^',ecolor='r', markeredgecolor='red', capthick=0.5, label=r"ALICE for $M_{X} < 200$ GeV")
p_isr = plt.errorbar(x8,y8,yerr=ery8,fmt='ko',ecolor='k', capthick=0.5, label=r"ISR for $M_{X}^{2}/s < 0.05$")
p_ua5 = plt.errorbar(x4,y4,yerr=ery4,fmt='w^',ecolor='b', markeredgecolor='blue', capthick=0.5, label=r"UA5 for $M_{X}^{2}/s < 0.05$")
p_cdf = plt.errorbar(x5,y5,yerr=ery5,fmt='wv',ecolor='b', markeredgecolor='blue', capthick=0.5, label=r"CDF for $M_{X}^{2}/s < 0.2$")
p_e710 = plt.errorbar(x6,y6,yerr=ery6,fmt='wo',ecolor='b', markeredgecolor='blue', capthick=0.5, label=r"E710 for $M_{X}^{2}/s < 0.05$")
p_ua4 = plt.errorbar(x7,y7,yerr=ery7,fmt='ws',ecolor='b', markeredgecolor='blue', capthick=0.5, label=r"UA4 for $M_{X}^{2}/s < 0.05$")


plt.plot(px1,py1, 'k-', linewidth=1.25)
plt.plot(px2,py2, 'k-', linewidth=1.25)
plt.plot(px3,py3, 'r--', linewidth=1.25)
plt.plot(px4,py4, 'r--', linewidth=1.25)
plt.plot(px5,py5, 'b:', linewidth=1.25)
plt.plot(px6,py6, 'b:', linewidth=1.25)

graph_data = [p_totem,p_cms,p_alice,p_isr,p_ua5,p_cdf,p_e710,p_ua4]
graph_data_label = ["TOTEM for 3.4 GeV $< M_{SD} <$ 1100 GeV",
  "CMS for $M_{X}^{2}/s < 0.05$ or $M_{Y}^{2}/s < 0.05$", "ALICE for $M_{X} < 200$ GeV",
  "ISR for $M_{X}^{2}/s < 0.05$",
  "UA5 for $M_{X}^{2}/s < 0.05$", "CDF for $M_{X}^{2}/s < 0.2$",
  "E710 for $M_{X}^{2}/s < 0.05$", "UA4 for $M_{X}^{2}/s < 0.05$"]

handle_curve_CT14 = mlines.Line2D([],[],color='black',linestyle='solid', linewidth=1.25, label=r'CT14')
handle_curve_CTE6L = mlines.Line2D([],[],color='red',linestyle='dashed', linewidth=1.25, label=r'NNPDF31')
handle_curve_MMHT = mlines.Line2D([],[],color='blue',linestyle='dotted', linewidth=1.25, label=r'MMHT')

artists=[handle_curve_CT14,handle_curve_CTE6L,handle_curve_MMHT]

leg1 = plt.legend(graph_data,graph_data_label,loc=2, ncol=1, shadow=False, fancybox=False, frameon=False, numpoints=1)
leg2 = plt.legend(handles=artists,loc=1, ncol=1, shadow=False, fancybox=False, frameon=False, numpoints=1) # this removes leg1 from the axes.
plt.gca().add_artist(leg1) # add leg1 as a separate artist to the axes

plt.fill_between(px3,py31,py32, color='r',alpha=0.2)
plt.fill_between(px4,py41,py42, color='r',alpha=0.2)
plt.fill_between(px5,py51,py52, color='b',alpha=0.2)
plt.fill_between(px6,py61,py62, color='b',alpha=0.2)

plt.xlabel(r"$\sqrt{s}$ [GeV]", fontdict=font)
plt.ylabel(r"$\sigma_{diff}(s)$ [mb]", fontdict=font)

plt.savefig("XSecDiff_BSG20orig_w002bPDFs_v3band.pdf")
#plt.show()

