import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.ndimage.filters import gaussian_filter
from matplotlib.cm import ScalarMappable


sims_per_point = 50
cm_k = plt.cm.get_cmap('RdBu')
cm_Ps = plt.cm.get_cmap('gist_heat')

#List fstem values run (set up to plot three diff values across panels)
f_stem_1 = 0.25
f_stem_2 = 0.5
f_stem_3 = 0.75

#Set colorbar attributes for eta and beta
eta_min=(10.0 - 50.0)
eta_max=(90.0 - 50.0)
color_ticks_eta = [-40.0 , -30.0, -20.0, -10.0, 0.0, 10.0, 20.0, 30.0, 40.0]
levels_eta=20
level_boundaries_eta = np.linspace(eta_min, eta_max, levels_eta + 1, endpoint=True)
sigma_eta = 5.0 #for Gaussian filter
stem_eta_contours = [-30, -20, -10, 0, 10, 20, 30]

beta_min=-.55
beta_max=0.01
color_ticks_beta = [-.5, -.4, -.3, -.2, -.1, 0.0]
levels_beta=20
level_boundaries_beta = np.linspace(beta_min, beta_max, levels_beta + 1, endpoint=True)
sigma_beta = 5.0 #for Gaussian filter
total_beta_contours = [-0.4, -0.3, -0.2, -0.1]


###Plot % conductance when Ps = mean P50 [eta]

filename_1 = "Output/varying_segmentation_n_{0}_fstem_{1}".format(int(sims_per_point), int(f_stem_1*100.0))
stem_P50_list_1 = np.loadtxt('%s/stem_P50_vals.txt' % (filename_1))
leaf_P50_list_1 = np.loadtxt('%s/leaf_P50_vals.txt' % (filename_1))
mean_k_frac_1 = np.loadtxt('%s/k_plant_frac_mean.txt' % (filename_1))
mean_stem_k_frac_1 = np.loadtxt('%s/k_stem_frac_mean.txt' % (filename_1))

filename_2 = "Output/varying_segmentation_n_{0}_fstem_{1}".format(int(sims_per_point), int(f_stem_2*100.0))
stem_P50_list_2 = np.loadtxt('%s/stem_P50_vals.txt' % (filename_2))
leaf_P50_list_2 = np.loadtxt('%s/leaf_P50_vals.txt' % (filename_2))
mean_k_frac_2 = np.loadtxt('%s/k_plant_frac_mean.txt' % (filename_2))
mean_stem_k_frac_2 = np.loadtxt('%s/k_stem_frac_mean.txt' % (filename_2))

filename_3 = "Output/varying_segmentation_n_{0}_fstem_{1}".format(int(sims_per_point), int(f_stem_3*100.0))
stem_P50_list_3 = np.loadtxt('%s/stem_P50_vals.txt' % (filename_3))
leaf_P50_list_3 = np.loadtxt('%s/leaf_P50_vals.txt' % (filename_3))
mean_k_frac_3 = np.loadtxt('%s/k_plant_frac_mean.txt' % (filename_3))
mean_stem_k_frac_3 = np.loadtxt('%s/k_stem_frac_mean.txt' % (filename_3))

meanP50_list_1 = [(a + b) / 2.0 for a,b in zip(stem_P50_list_1, leaf_P50_list_1)]
meanP50_list_1 = np.around(meanP50_list_1, 3)
meanP50_list_2 = [(a + b) / 2.0 for a,b in zip(stem_P50_list_2, leaf_P50_list_2)]
meanP50_list_2 = np.around(meanP50_list_2, 3)
meanP50_list_3 = [(a + b) / 2.0 for a,b in zip(stem_P50_list_3, leaf_P50_list_3)]
meanP50_list_3 = np.around(meanP50_list_3, 3)
dp50_list_1 = np.around(leaf_P50_list_1, 4) - np.around(stem_P50_list_1, 4)
dp50_list_1 = np.around(dp50_list_1, 2)
dp50_list_2 = np.around(leaf_P50_list_2, 4) - np.around(stem_P50_list_2, 4)
dp50_list_2 = np.around(dp50_list_2, 2)
dp50_list_3 = np.around(leaf_P50_list_3, 4) - np.around(stem_P50_list_3, 4)
dp50_list_3 = np.around(dp50_list_3, 2)


##Gaussian Filtering for k_stem
mean_stem_k_frac_1 = (mean_stem_k_frac_1*100.0)-50.0
df_1= pd.DataFrame(dict(x1=meanP50_list_1, y1=dp50_list_1, z1=mean_stem_k_frac_1))
xcol_1, ycol_1, zcol_1 = "x1", "y1", "z1"
df_1 = df_1.sort_values(by=[xcol_1, ycol_1])
xvals_1 = df_1[xcol_1].unique()
yvals_1 = df_1[ycol_1].unique()
zvals_1 = df_1[zcol_1].values.reshape(len(xvals_1), len(yvals_1)).T
data2_1 = gaussian_filter(zvals_1, sigma_eta)

mean_stem_k_frac_2 = (mean_stem_k_frac_2*100.0)-50.0
df_2= pd.DataFrame(dict(x2=meanP50_list_2, y2=dp50_list_2, z2=mean_stem_k_frac_2))
xcol_2, ycol_2, zcol_2 = "x2", "y2", "z2"
df_2 = df_2.sort_values(by=[xcol_2, ycol_2])
xvals_2 = df_2[xcol_2].unique()
yvals_2 = df_2[ycol_2].unique()
zvals_2 = df_2[zcol_2].values.reshape(len(xvals_2), len(yvals_2)).T
data2_2 = gaussian_filter(zvals_2, sigma_eta)

mean_stem_k_frac_3 = (mean_stem_k_frac_3*100.0)-50.0
df_3= pd.DataFrame(dict(x3=meanP50_list_3, y3=dp50_list_3, z3=mean_stem_k_frac_3))
xcol_3, ycol_3, zcol_3 = "x3", "y3", "z3"
df_3 = df_3.sort_values(by=[xcol_3, ycol_3])
xvals_3 = df_3[xcol_3].unique()
yvals_3 = df_3[ycol_3].unique()
zvals_3 = df_3[zcol_3].values.reshape(len(xvals_3), len(yvals_3)).T
data2_3 = gaussian_filter(zvals_3, sigma_eta)


##Plot eta for stem
fig, axs = plt.subplots(1, 3, figsize=(9,4), sharey=True)

CS1= axs[0].contourf(xvals_1, yvals_1, data2_1, levels_eta, cmap=cm_k, vmin=eta_min, vmax=eta_max)
CS1_1=axs[0].contour(xvals_1, yvals_1, data2_1, levels=stem_eta_contours, colors='k')
axs[0].clabel(CS1_1, inline=1, fontsize=10, fmt='%.f')
axs[0].set_ylim(min(yvals_1), max(yvals_1))
axs[0].set_title(r'$f_{stem,0}$ = 0.25')

CS2= axs[1].contourf(xvals_2, yvals_2, data2_2, levels_eta, cmap=cm_k, vmin=eta_min, vmax=eta_max)
CS2_1=axs[1].contour(xvals_2, yvals_2, data2_2, levels=stem_eta_contours, colors='k')
axs[1].clabel(CS2_1, inline=1, fontsize=10, fmt='%.f')
axs[1].set_ylim(min(yvals_2), max(yvals_2))
axs[1].set_title(r'$f_{stem,0}$ = 0.50')


CS3= axs[2].contourf(xvals_3, yvals_3, data2_3, levels_eta, cmap=cm_k, vmin=eta_min, vmax=eta_max)
CS3_1=axs[2].contour(xvals_3, yvals_3, data2_3, levels=stem_eta_contours, colors='k')
axs[2].clabel(CS3_1, inline=1, fontsize=10, fmt='%.f')
axs[2].set_ylim(min(yvals_3), max(yvals_3))
axs[2].set_title(r'$f_{stem,0}$ = 0.75')
CS3_map = ScalarMappable(norm=CS3.norm, cmap=CS3.cmap)

for ax in axs.flat:
    ax.set(xlabel=r'Mean $P_{50}$ (MPa)', ylabel=r'$\Delta P_{50}$ (MPa)')

for ax in axs.flat:
    ax.label_outer()

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.70])
CS3_map.set_array([])
fig.colorbar(
    CS3_map,
    boundaries=level_boundaries_eta,
    values=(level_boundaries_eta[:-1] + level_boundaries_eta[1:]) / 2,
    ticks=color_ticks_eta,
    cax=cbar_ax,
    label=r'$\eta$ (%)'
)

plt.savefig('Other_Plots/experiment_4_eta_stem.pdf')
plt.show()


###Plot Psoil value when k=50% [Beta]

filename_1 = "Output/varying_segmentation_n_{0}_fstem_{1}".format(int(sims_per_point), int(f_stem_1*100.0))
stem_P50_list_1 = np.loadtxt('%s/stem_P50_vals.txt' % (filename_1))
leaf_P50_list_1 = np.loadtxt('%s/leaf_P50_vals.txt' % (filename_1))
mean_total_Ps_1 = np.loadtxt('%s/total50_Ps.txt' % (filename_1))
mean_stem_Ps_1 = np.loadtxt('%s/stem50_Ps.txt' % (filename_1))
mean_P50_list_1 = [(a + b) / 2.0 for a,b in zip(stem_P50_list_1, leaf_P50_list_1)]
mean_P50_list_1 = np.around(mean_P50_list_1, 3)
Ps_margin_1 = mean_P50_list_1 - mean_total_Ps_1
Ps_stem_margin_1 = stem_P50_list_1 - mean_stem_Ps_1

filename_2 = "Output/varying_segmentation_n_{0}_fstem_{1}".format(int(sims_per_point), int(f_stem_2*100.0))
stem_P50_list_2 = np.loadtxt('%s/stem_P50_vals.txt' % (filename_2))
leaf_P50_list_2 = np.loadtxt('%s/leaf_P50_vals.txt' % (filename_2))
mean_total_Ps_2 = np.loadtxt('%s/total50_Ps.txt' % (filename_2))
mean_stem_Ps_2 = np.loadtxt('%s/stem50_Ps.txt' % (filename_2))
mean_P50_list_2 = [(a + b) / 2.0 for a,b in zip(stem_P50_list_2, leaf_P50_list_2)]
mean_P50_list_2 = np.around(mean_P50_list_2, 3)
Ps_margin_2 = mean_P50_list_2 - mean_total_Ps_2
Ps_stem_margin_2 = stem_P50_list_2 - mean_stem_Ps_2

filename_3 = "Output/varying_segmentation_n_{0}_fstem_{1}".format(int(sims_per_point), int(f_stem_3*100.0))
stem_P50_list_3 = np.loadtxt('%s/stem_P50_vals.txt' % (filename_3))
leaf_P50_list_3 = np.loadtxt('%s/leaf_P50_vals.txt' % (filename_3))
mean_total_Ps_3 = np.loadtxt('%s/total50_Ps.txt' % (filename_3))
mean_stem_Ps_3 = np.loadtxt('%s/stem50_Ps.txt' % (filename_3))
mean_P50_list_3 = [(a + b) / 2.0 for a,b in zip(stem_P50_list_3, leaf_P50_list_3)]
mean_P50_list_3 = np.around(mean_P50_list_3, 3)
Ps_margin_3 = mean_P50_list_3 - mean_total_Ps_3
Ps_stem_margin_3 = stem_P50_list_3 - mean_stem_Ps_3

dp50_list_1 = np.around(leaf_P50_list_1, 4) - np.around(stem_P50_list_1, 4)
dp50_list_1 = np.around(dp50_list_1, 2)
dp50_list_2 = np.around(leaf_P50_list_2, 4) - np.around(stem_P50_list_2, 4)
dp50_list_2 = np.around(dp50_list_2, 2)
dp50_list_3 = np.around(leaf_P50_list_3, 4) - np.around(stem_P50_list_3, 4)
dp50_list_3 = np.around(dp50_list_3, 2)

##Gaussian Filtering

#first f_stem value

df_1= pd.DataFrame(dict(x1=mean_P50_list_1, y1=dp50_list_1, z1=Ps_margin_1))
xcol_1, ycol_1, zcol_1 = "x1", "y1", "z1"
df_1 = df_1.sort_values(by=[xcol_1, ycol_1])
xvals_1 = df_1[xcol_1].unique()
yvals_1 = df_1[ycol_1].unique()
zvals_1 = df_1[zcol_1].values.reshape(len(xvals_1), len(yvals_1)).T
data2_1 = gaussian_filter(zvals_1, sigma_beta)

#2nd f_stem value

df_2= pd.DataFrame(dict(x2=mean_P50_list_2, y2=dp50_list_2, z2=Ps_margin_2))
xcol_2, ycol_2, zcol_2 = "x2", "y2", "z2"
df_2 = df_2.sort_values(by=[xcol_2, ycol_2])
xvals_2 = df_2[xcol_2].unique()
yvals_2 = df_2[ycol_2].unique()
zvals_2 = df_2[zcol_2].values.reshape(len(xvals_2), len(yvals_2)).T
data2_2 = gaussian_filter(zvals_2, sigma_beta)

#3rd f_stem value

df_3= pd.DataFrame(dict(x3=mean_P50_list_3, y3=dp50_list_3, z3=Ps_margin_3))
xcol_3, ycol_3, zcol_3 = "x3", "y3", "z3"
df_3 = df_3.sort_values(by=[xcol_3, ycol_3])
xvals_3 = df_3[xcol_3].unique()
yvals_3 = df_3[ycol_3].unique()
zvals_3 = df_3[zcol_3].values.reshape(len(xvals_3), len(yvals_3)).T
data2_3 = gaussian_filter(zvals_3, sigma_beta)

##Plot
fig, axs = plt.subplots(1, 3, figsize=(9,4), sharey=True)

CS1= axs[0].contourf(xvals_1, yvals_1, data2_1, levels_beta, cmap=cm_Ps, vmin=beta_min, vmax=beta_max)
CS1_1=axs[0].contour(xvals_1, yvals_1, data2_1, levels=total_beta_contours, colors='k')
axs[0].clabel(CS1_1, inline=1, fontsize=10, fmt='%.1f')
axs[0].set_ylim(min(yvals_1), max(yvals_1))
axs[0].set_title(r'$f_{stem,0}$ = 0.25')

CS2= axs[1].contourf(xvals_2, yvals_2, data2_2, levels_beta, cmap=cm_Ps, vmin=beta_min, vmax=beta_max)
CS2_1=axs[1].contour(xvals_2, yvals_2, data2_2, levels=total_beta_contours, colors='k')
axs[1].clabel(CS2_1, inline=1, fontsize=10, fmt='%.1f')
axs[1].set_ylim(min(yvals_2), max(yvals_2))
axs[1].set_title(r'$f_{stem,0}$ = 0.50')


CS3= axs[2].contourf(xvals_3, yvals_3, data2_3, levels_beta, cmap=cm_Ps, vmin=beta_min, vmax=beta_max)
CS3_1=axs[2].contour(xvals_3, yvals_3, data2_3, levels=total_beta_contours, colors='k')
axs[2].clabel(CS3_1, inline=1, fontsize=10, fmt='%.1f')
axs[2].set_ylim(min(yvals_3), max(yvals_3))
axs[2].set_title(r'$f_{stem,0}$ = 0.75')
CS3_map = ScalarMappable(norm=CS3.norm, cmap=CS3.cmap)

for ax in axs.flat:
    ax.set(xlabel=r'Mean $P_{50}$ (MPa)', ylabel=r'$\Delta P_{50}$ (MPa)')

for ax in axs.flat:
    ax.label_outer()

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.70])
CS3_map.set_array([])
fig.colorbar(
    CS3_map,
    boundaries=level_boundaries_beta,
    values=(level_boundaries_beta[:-1] + level_boundaries_beta[1:]) / 2,
    ticks=color_ticks_beta,
    cax=cbar_ax,
    label=r'$\beta$ (MPa)'
)

plt.savefig('Other_Plots/experiment_4_beta_total.pdf')
plt.show()


