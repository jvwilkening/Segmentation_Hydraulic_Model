import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.ndimage.filters import gaussian_filter
from matplotlib.cm import ScalarMappable

sims_per_point = 50
gmin_frac = .05
cm_diff = plt.cm.get_cmap('PiYG')

#List magnitude of delta P50 values run (set up to plot three diff values across panels)
dp50_1 = 0.1
dp50_2 = 0.5
dp50_3 = 1.0

#max and min mean P50 range
max_mean_P50 = -2.0
min_mean_P50 = -6.0

#set range of lambda total and colorbar attributes
total_min_diff=-55.0
total_max_diff=55.0
color_ticks_total = [-50.0, -25.0, 0.0, 25.0, 50.0]
total_diff_contours = [-30.0, -20.0, -10.0, 0.0, 10.0, 20.0, 30.0, 40.0]
levels_total=20
level_boundaries_total = np.linspace(total_min_diff, total_max_diff, levels_total+1, endpoint=True)
sigma_total = 5.0 #used for Gaussian filter


#set range of lambda stem and colorbar attributes
stem_min_diff=-65.0
stem_max_diff=65.0
color_ticks_stem = [-50.0, -25.0, 0.0, 25.0, 50.0]
stem_diff_contours = [0.0, 25.0, 50.0]
#color_ticks_total = np.linspace(total_min_k, total_max_k, 5)
levels_stem=20
level_boundaries_stem = np.linspace(total_min_diff, total_max_diff, levels_stem+1, endpoint=True)
sigma_stem = 5.0 #used for Gaussian filter

###Plot % diff between vulnerability segmentation patterns when Ps = max(P50) in plant (ie most vulnerable)
filename_1 = "Output/P50_vs_fstem_n_{0}_dp50_{2}".format(int(sims_per_point), int(dp50_1*100.0))
rev_stem_1 = np.loadtxt('%s/rev_k_stem_frac_mean.txt' % (filename_1))
conv_stem_1 = np.loadtxt('%s/conv_k_stem_frac_mean.txt' % (filename_1))
conv_stem_ci_1 = np.loadtxt('%s/conv_k_stem_frac_ci.txt' % (filename_1))
rev_stem_ci_1 = np.loadtxt('%s/rev_k_stem_frac_ci.txt' % (filename_1))
rev_plant_1 = np.loadtxt('%s/rev_k_plant_frac_mean.txt' % (filename_1))
conv_plant_1 = np.loadtxt('%s/conv_k_plant_frac_mean.txt' % (filename_1))
conv_plant_ci_1 = np.loadtxt('%s/conv_k_plant_frac_ci.txt' % (filename_1))
rev_plant_ci_1 = np.loadtxt('%s/rev_k_plant_frac_ci.txt' % (filename_1))
max_P50_vals_1 = np.loadtxt('%s/max_P50_vals.txt' % (filename_1))
f_stem_vals_1 = np.loadtxt('%s/f_stem_vals.txt' % (filename_1))
lambda_total_1 = (conv_plant_1 - rev_plant_1) * 100.0
lambda_stem_1 = (conv_stem_1 - rev_stem_1) * 100.0
mean_P50_vals_1 = (max_P50_vals_1 + (max_P50_vals_1 - dp50_1))/2.0

filename_2 = "Output/P50_vs_fstem_n_{0}_gfrac_{1}_dp50_{2}".format(int(sims_per_point), int(gmin_frac*100.0), int(dp50_2*100.0))
rev_stem_2 = np.loadtxt('%s/rev_k_stem_frac_mean.txt' % (filename_2))
conv_stem_2 = np.loadtxt('%s/conv_k_stem_frac_mean.txt' % (filename_2))
conv_stem_ci_2 = np.loadtxt('%s/conv_k_stem_frac_ci.txt' % (filename_2))
rev_stem_ci_2 = np.loadtxt('%s/rev_k_stem_frac_ci.txt' % (filename_2))
rev_plant_2 = np.loadtxt('%s/rev_k_plant_frac_mean.txt' % (filename_2))
conv_plant_2 = np.loadtxt('%s/conv_k_plant_frac_mean.txt' % (filename_2))
conv_plant_ci_2 = np.loadtxt('%s/conv_k_plant_frac_ci.txt' % (filename_2))
rev_plant_ci_2 = np.loadtxt('%s/rev_k_plant_frac_ci.txt' % (filename_2))
max_P50_vals_2 = np.loadtxt('%s/max_P50_vals.txt' % (filename_2))
f_stem_vals_2 = np.loadtxt('%s/f_stem_vals.txt' % (filename_2))
lambda_total_2 = (conv_plant_2 - rev_plant_2) * 100.0
lambda_stem_2 = (conv_stem_2 - rev_stem_2) * 100.0
mean_P50_vals_2 = (max_P50_vals_2 + (max_P50_vals_2 - dp50_2))/2.0

filename_3 = "Output/P50_vs_fstem_n_{0}_gfrac_{1}_dp50_{2}".format(int(sims_per_point), int(gmin_frac*100.0), int(dp50_3*100.0))
rev_stem_3 = np.loadtxt('%s/rev_k_stem_frac_mean.txt' % (filename_3))
conv_stem_3 = np.loadtxt('%s/conv_k_stem_frac_mean.txt' % (filename_3))
conv_stem_ci_3 = np.loadtxt('%s/conv_k_stem_frac_ci.txt' % (filename_3))
rev_stem_ci_3 = np.loadtxt('%s/rev_k_stem_frac_ci.txt' % (filename_3))
rev_plant_3 = np.loadtxt('%s/rev_k_plant_frac_mean.txt' % (filename_3))
conv_plant_3 = np.loadtxt('%s/conv_k_plant_frac_mean.txt' % (filename_3))
conv_plant_ci_3 = np.loadtxt('%s/conv_k_plant_frac_ci.txt' % (filename_3))
rev_plant_ci_3 = np.loadtxt('%s/rev_k_plant_frac_ci.txt' % (filename_3))
max_P50_vals_3 = np.loadtxt('%s/max_P50_vals.txt' % (filename_3))
f_stem_vals_3 = np.loadtxt('%s/f_stem_vals.txt' % (filename_3))
lambda_total_3 = (conv_plant_3 - rev_plant_3) * 100.0
lambda_stem_3 = (conv_stem_3 - rev_stem_3) * 100.0
mean_P50_vals_3 = (max_P50_vals_3 + (max_P50_vals_3 - dp50_3))/2.0

##Gaussian Filtering

df_1= pd.DataFrame(dict(x1=mean_P50_vals_1, y1=f_stem_vals_1, z1=lambda_total_1))
xcol_1, ycol_1, zcol_1 = "x1", "y1", "z1"
df_1 = df_1.sort_values(by=[xcol_1, ycol_1])
xvals_1 = df_1[xcol_1].unique()
yvals_1 = df_1[ycol_1].unique()
zvals_1 = df_1[zcol_1].values.reshape(len(xvals_1), len(yvals_1)).T
data2_1 = gaussian_filter(zvals_1, sigma_total)

df_2= pd.DataFrame(dict(x2=mean_P50_vals_2, y2=f_stem_vals_2, z2=lambda_total_2))
xcol_2, ycol_2, zcol_2 = "x2", "y2", "z2"
df_2 = df_2.sort_values(by=[xcol_2, ycol_2])
xvals_2 = df_2[xcol_2].unique()
yvals_2 = df_2[ycol_2].unique()
zvals_2 = df_2[zcol_2].values.reshape(len(xvals_2), len(yvals_2)).T
data2_2 = gaussian_filter(zvals_2, sigma_total)

df_3= pd.DataFrame(dict(x3=mean_P50_vals_3, y3=f_stem_vals_3, z3=lambda_total_3))
xcol_3, ycol_3, zcol_3 = "x3", "y3", "z3"
df_3 = df_3.sort_values(by=[xcol_3, ycol_3])
xvals_3 = df_3[xcol_3].unique()
yvals_3 = df_3[ycol_3].unique()
zvals_3 = df_3[zcol_3].values.reshape(len(xvals_3), len(yvals_3)).T
data2_3 = gaussian_filter(zvals_3, sigma_total)

#Plot data
fig, axs = plt.subplots(1, 3, figsize=(9,4), sharey=True)

CS1= axs[0].contourf(xvals_1, yvals_1, data2_1, levels_total, cmap=cm_diff, vmin=total_min_diff, vmax=total_max_diff)
CS1_1=axs[0].contour(xvals_1, yvals_1, data2_1, levels=total_diff_contours, colors='k')
axs[0].clabel(CS1_1, inline=1, fontsize=10, fmt='%.0f')
axs[0].hlines(0.0, min(xvals_1), max(xvals_1), linestyles='dashed')
axs[0].set_ylim(min(yvals_1), max(yvals_1))
axs[0].set_xlim(min_mean_P50, max_mean_P50)
axs[0].set_title(r'$|\Delta P_{50}|$ = 0.1 MPa')

CS2= axs[1].contourf(xvals_2, yvals_2, data2_2, levels_total, cmap=cm_diff, vmin=total_min_diff, vmax=total_max_diff)
CS2_1=axs[1].contour(xvals_2, yvals_2, data2_2, levels=total_diff_contours, colors='k')
axs[1].clabel(CS2_1, inline=1, fontsize=10, fmt='%.0f')
axs[1].hlines(0.0, min(xvals_2), max(xvals_2), linestyles='dashed')
axs[1].set_ylim(min(yvals_2), max(yvals_2))
axs[1].set_xlim(min_mean_P50, max_mean_P50)
axs[1].set_title(r'$|\Delta P_{50}|$ = 0.5 MPa')

CS3= axs[2].contourf(xvals_3, yvals_3, data2_3, levels_total, cmap=cm_diff, vmin=total_min_diff, vmax=total_max_diff)
CS3_1=axs[2].contour(xvals_3, yvals_3, data2_3, levels=total_diff_contours, colors='k')
axs[2].clabel(CS3_1, inline=1, fontsize=10, fmt='%.0f')
axs[2].hlines(0.0, min(xvals_3), max(xvals_3), linestyles='dashed')
axs[2].set_ylim(min(yvals_3), max(yvals_3))
axs[2].set_xlim(min_mean_P50, max_mean_P50)
axs[2].set_title(r'$|\Delta P_{50}|$ = 1.0 MPa')

CS3_map = ScalarMappable(norm=CS3.norm, cmap=CS3.cmap)

for ax in axs.flat:
    ax.set(xlabel=r'Mean $P_{50}$ in Plant (MPa)', ylabel=r'$f_{stem,0}$')

for ax in axs.flat:
    ax.label_outer()

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.70])
CS3_map.set_array([])
fig.colorbar(
    CS3_map,
    boundaries=level_boundaries_total,
    values=(level_boundaries_total[:-1] + level_boundaries_total[1:]) / 2,
    ticks=color_ticks_total,
    cax=cbar_ax,
    label=r'$\lambda_{total}$ (%)'
)

plt.savefig('Other_Plots/experiment_2_lambda_total.pdf')
plt.show()


##Plot lambda stem results

#Gaussian Filtering

df_1= pd.DataFrame(dict(x1=mean_P50_vals_1, y1=f_stem_vals_1, z1=lambda_stem_1))
xcol_1, ycol_1, zcol_1 = "x1", "y1", "z1"
df_1 = df_1.sort_values(by=[xcol_1, ycol_1])
xvals_1 = df_1[xcol_1].unique()
yvals_1 = df_1[ycol_1].unique()
zvals_1 = df_1[zcol_1].values.reshape(len(xvals_1), len(yvals_1)).T
data2_1 = gaussian_filter(zvals_1, sigma_stem)

df_2= pd.DataFrame(dict(x2=mean_P50_vals_2, y2=f_stem_vals_2, z2=lambda_stem_2))
xcol_2, ycol_2, zcol_2 = "x2", "y2", "z2"
df_2 = df_2.sort_values(by=[xcol_2, ycol_2])
xvals_2 = df_2[xcol_2].unique()
yvals_2 = df_2[ycol_2].unique()
zvals_2 = df_2[zcol_2].values.reshape(len(xvals_2), len(yvals_2)).T
data2_2 = gaussian_filter(zvals_2, sigma_stem)

df_3= pd.DataFrame(dict(x3=mean_P50_vals_3, y3=f_stem_vals_3, z3=lambda_stem_3))
xcol_3, ycol_3, zcol_3 = "x3", "y3", "z3"
df_3 = df_3.sort_values(by=[xcol_3, ycol_3])
xvals_3 = df_3[xcol_3].unique()
yvals_3 = df_3[ycol_3].unique()
zvals_3 = df_3[zcol_3].values.reshape(len(xvals_3), len(yvals_3)).T
data2_3 = gaussian_filter(zvals_3, sigma_stem)

#Plot data
fig, axs = plt.subplots(1, 3, figsize=(9,4), sharey=True)

CS1= axs[0].contourf(xvals_1, yvals_1, data2_1, levels_stem, cmap=cm_diff, vmin=stem_min_diff, vmax=stem_max_diff)
CS1_1=axs[0].contour(xvals_1, yvals_1, data2_1, levels=stem_diff_contours, colors='k')
axs[0].clabel(CS1_1, inline=1, fontsize=10, fmt='%.0f')
axs[0].hlines(0.0, min(xvals_1), max(xvals_1), linestyles='dashed')
axs[0].set_ylim(min(yvals_1), max(yvals_1))
axs[0].set_xlim(min_mean_P50, max_mean_P50)
axs[0].set_title(r'$|\Delta P_{50}|$ = 0.1 MPa')

CS2= axs[1].contourf(xvals_2, yvals_2, data2_2, levels_stem, cmap=cm_diff, vmin=stem_min_diff, vmax=stem_max_diff)
CS2_1=axs[1].contour(xvals_2, yvals_2, data2_2, levels=stem_diff_contours, colors='k')
axs[1].clabel(CS2_1, inline=1, fontsize=10, fmt='%.0f')
axs[1].hlines(0.0, min(xvals_2), max(xvals_2), linestyles='dashed')
axs[1].set_ylim(min(yvals_2), max(yvals_2))
axs[1].set_xlim(min_mean_P50, max_mean_P50)
axs[1].set_title(r'$|\Delta P_{50}|$ = 0.5 MPa')

CS3= axs[2].contourf(xvals_3, yvals_3, data2_3, levels_stem, cmap=cm_diff, vmin=stem_min_diff, vmax=stem_max_diff)
#CS3_1=axs[2].contour(xvals_3, yvals_3, data2_3, levels=stem_diff_contours, colors='k')
#axs[2].clabel(CS3_1, inline=1, fontsize=10, fmt='%.0f')
axs[2].hlines(0.0, min(xvals_3), max(xvals_3), linestyles='dashed')
axs[2].set_ylim(min(yvals_3), max(yvals_3))
axs[2].set_xlim(min_mean_P50, max_mean_P50)
axs[2].set_title(r'$|\Delta P_{50}|$ = 1.0 MPa')

CS3_map = ScalarMappable(norm=CS3.norm, cmap=CS3.cmap)

for ax in axs.flat:
    ax.set(xlabel=r'Mean $P_{50}$ in Plant (MPa)', ylabel=r'$f_{stem,0}$')

for ax in axs.flat:
    ax.label_outer()

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.70])
CS3_map.set_array([])
fig.colorbar(
    CS3_map,
    boundaries=level_boundaries_stem,
    values=(level_boundaries_stem[:-1] + level_boundaries_stem[1:]) / 2,
    ticks=color_ticks_stem,
    cax=cbar_ax,
    label=r'$\lambda_{stem}$ (%)'
)

plt.savefig('Other_Plots/experiment_2_lambda_stem.pdf')
plt.show()