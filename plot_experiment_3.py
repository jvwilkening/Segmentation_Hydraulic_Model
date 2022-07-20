import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.ndimage.filters import gaussian_filter
from matplotlib.cm import ScalarMappable

sims_per_point = 50
cm_Ps = plt.cm.get_cmap('gist_heat')
cm_k = plt.cm.get_cmap('RdBu')

#List mean P50 values run (set up to plot three diff values across panels)
P50_mean_1 = 2.0
P50_mean_2 = 4.0
P50_mean_3 = 6.0


#Sets colorbar attributes for eta and beta
eta_min= 5.0 - 50.0
eta_max= 90.0 - 50.0
color_ticks_eta = [-40.0, -30.0, -20.0, -10.0, 0.0, 10.0, 20.0, 30.0, 40.0]
eta_contours = [-25.0, 00.0, 25.0]
levels_eta=20
level_boundaries_eta = np.linspace(eta_min, eta_max, levels_eta + 1, endpoint=True)
sigma_eta = 5.0 #used for Gaussian filter

beta_min=-.7
beta_max=0.15
color_ticks_beta = [-.7, -.6, -.5, -.4, -.3, -.2, -.1, 0.0, 0.1]
beta_contours = [-.6, -.4, -.2, 0.0]
levels_beta=20
level_boundaries_beta = np.linspace(beta_min, beta_max, levels_beta + 1, endpoint=True)
sigma_beta = 5.0 #used for Gaussian filter



###Plot % total conductance when Ps = mean P50 [eta]
filename_1 = "Output/varying_segmentation_and_fstem_n_{0}_P50mean_{1}".format(int(sims_per_point), int(P50_mean_1))
f_stem_vals_1 = np.loadtxt('%s/fstem_vals.txt' % (filename_1))
dP50_vals_1 = np.loadtxt('%s/dP50_vals.txt' % (filename_1))
mean_k_frac_1 = np.loadtxt('%s/k_plant_frac_mean.txt' % (filename_1))
eta_1 = np.loadtxt('%s/k_stem_frac_mean.txt' % (filename_1))
eta_1 = eta_1 * 100.0 - 50.0


filename_2 = "Output/varying_segmentation_and_fstem_n_{0}_P50mean_{1}".format(int(sims_per_point), int(P50_mean_2))
f_stem_vals_2 = np.loadtxt('%s/fstem_vals.txt' % (filename_2))
dP50_vals_2 = np.loadtxt('%s/dP50_vals.txt' % (filename_2))
mean_k_frac_2 = np.loadtxt('%s/k_plant_frac_mean.txt' % (filename_2))
eta_2 = np.loadtxt('%s/k_stem_frac_mean.txt' % (filename_2))
eta_2 = eta_2 * 100.0 - 50.0

filename_3 = "Output/varying_segmentation_and_fstem_n_{0}_P50mean_{1}".format(int(sims_per_point), int(P50_mean_3))
f_stem_vals_3 = np.loadtxt('%s/fstem_vals.txt' % (filename_3))
dP50_vals_3 = np.loadtxt('%s/dP50_vals.txt' % (filename_3))
mean_k_frac_3 = np.loadtxt('%s/k_plant_frac_mean.txt' % (filename_3))
eta_3 = np.loadtxt('%s/k_stem_frac_mean.txt' % (filename_3))
eta_3 = eta_3 * 100.0 - 50.0


##Plot eta for stem conductance results
##Gaussian Filtering

df_1= pd.DataFrame(dict(x1=dP50_vals_1, y1=f_stem_vals_1, z1=eta_1))
xcol_1, ycol_1, zcol_1 = "x1", "y1", "z1"
df_1 = df_1.sort_values(by=[xcol_1, ycol_1])
xvals_1 = df_1[xcol_1].unique()
yvals_1 = df_1[ycol_1].unique()
zvals_1 = df_1[zcol_1].values.reshape(len(xvals_1), len(yvals_1)).T
data2_1 = gaussian_filter(zvals_1, sigma_eta)

df_2= pd.DataFrame(dict(x2=dP50_vals_2, y2=f_stem_vals_2, z2=eta_2))
xcol_2, ycol_2, zcol_2 = "x2", "y2", "z2"
df_2 = df_2.sort_values(by=[xcol_2, ycol_2])
xvals_2 = df_2[xcol_2].unique()
yvals_2 = df_2[ycol_2].unique()
zvals_2 = df_2[zcol_2].values.reshape(len(xvals_2), len(yvals_2)).T
data2_2 = gaussian_filter(zvals_2, sigma_eta)

df_3= pd.DataFrame(dict(x3=dP50_vals_3, y3=f_stem_vals_3, z3=eta_3))
xcol_3, ycol_3, zcol_3 = "x3", "y3", "z3"
df_3 = df_3.sort_values(by=[xcol_3, ycol_3])
xvals_3 = df_3[xcol_3].unique()
yvals_3 = df_3[ycol_3].unique()
zvals_3 = df_3[zcol_3].values.reshape(len(xvals_3), len(yvals_3)).T
data2_3 = gaussian_filter(zvals_3, sigma_eta)

#Plot data
fig, axs = plt.subplots(1, 3, figsize=(9,4), sharey=True)

CS1= axs[0].contourf(xvals_1, yvals_1, data2_1, levels_eta, cmap=cm_k, vmin=eta_min, vmax=eta_max)
CS1_1=axs[0].contour(xvals_1, yvals_1, data2_1, levels=eta_contours, colors='k')
axs[0].clabel(CS1_1, inline=1, fontsize=10, fmt='%.0f')
axs[0].set_ylim(min(yvals_1), max(yvals_1))
axs[0].set_title(r'Mean $P_{50}$ = -2.0 MPa')

CS2= axs[1].contourf(xvals_2, yvals_2, data2_2, levels_eta, cmap=cm_k, vmin=eta_min, vmax=eta_max)
CS2_1=axs[1].contour(xvals_2, yvals_2, data2_2, levels=eta_contours, colors='k')
axs[1].clabel(CS2_1, inline=1, fontsize=10, fmt='%.0f')
axs[1].set_ylim(min(yvals_2), max(yvals_2))
axs[1].set_title(r'Mean $P_{50}$ = -4.0 MPa')

CS3= axs[2].contourf(xvals_3, yvals_3, data2_3, levels_eta, cmap=cm_k, vmin=eta_min, vmax=eta_max)
CS3_1=axs[2].contour(xvals_3, yvals_3, data2_3, levels=eta_contours, colors='k')
axs[2].clabel(CS3_1, inline=1, fontsize=10, fmt='%.0f')
axs[2].set_ylim(min(yvals_3), max(yvals_3))
axs[2].set_title(r'Mean $P_{50}$ = -6.0 MPa')

CS3_map = ScalarMappable(norm=CS3.norm, cmap=CS3.cmap)

for ax in axs.flat:
    ax.set(xlabel=r'$\Delta P_{50}$ (MPa)', ylabel=r'$f_{stem,0}$')

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
plt.savefig('Other_Plots/experiment_2_eta_stem.pdf')
plt.show()


###Plot beta for total conductance results
filename_1 = "Output/varying_segmentation_and_fstem_with_Ps_n_{0}_P50mean_{1}".format(int(sims_per_point), int(P50_mean_1))
f_stem_vals_1 = np.loadtxt('%s/fstem_vals.txt' % (filename_1))
dP50_vals_1 = np.loadtxt('%s/dP50_vals.txt' % (filename_1))
mean_total_Ps_1 = np.loadtxt('%s/total50_Ps.txt' % (filename_1))
mean_stem_Ps_1 = np.loadtxt('%s/stem50_Ps.txt' % (filename_1))
beta_1 = -P50_mean_1 - mean_total_Ps_1


filename_2 = "Output/varying_segmentation_and_fstem_with_Ps_n_{0}P50mean_{1}".format(int(sims_per_point),int(P50_mean_2))
f_stem_vals_2 = np.loadtxt('%s/fstem_vals.txt' % (filename_2))
dP50_vals_2 = np.loadtxt('%s/dP50_vals.txt' % (filename_2))
mean_total_Ps_2 = np.loadtxt('%s/total50_Ps.txt' % (filename_2))
mean_stem_Ps_2 = np.loadtxt('%s/stem50_Ps.txt' % (filename_2))
beta_2 = -P50_mean_2 - mean_total_Ps_2


filename_3 = "Output/varying_segmentation_and_fstem_with_Ps_n_{0}_P50mean_{1}_debug".format(int(sims_per_point), int(P50_mean_3))
f_stem_vals_3 = np.loadtxt('%s/fstem_vals.txt' % (filename_3))
dP50_vals_3 = np.loadtxt('%s/dP50_vals.txt' % (filename_3))
mean_total_Ps_3 = np.loadtxt('%s/total50_Ps.txt' % (filename_3))
mean_stem_Ps_3 = np.loadtxt('%s/stem50_Ps.txt' % (filename_3))
beta_3 = -P50_mean_3 - mean_total_Ps_3


#Gaussian Filtering
df_1= pd.DataFrame(dict(x1=dP50_vals_1, y1=f_stem_vals_1, z1=beta_1))
xcol_1, ycol_1, zcol_1 = "x1", "y1", "z1"
df_1 = df_1.sort_values(by=[xcol_1, ycol_1])
xvals_1 = df_1[xcol_1].unique()
yvals_1 = df_1[ycol_1].unique()
zvals_1 = df_1[zcol_1].values.reshape(len(xvals_1), len(yvals_1)).T
data2_1 = gaussian_filter(zvals_1, sigma_beta)

df_2= pd.DataFrame(dict(x2=dP50_vals_2, y2=f_stem_vals_2, z2=beta_2))
xcol_2, ycol_2, zcol_2 = "x2", "y2", "z2"
df_2 = df_2.sort_values(by=[xcol_2, ycol_2])
xvals_2 = df_2[xcol_2].unique()
yvals_2 = df_2[ycol_2].unique()
zvals_2 = df_2[zcol_2].values.reshape(len(xvals_2), len(yvals_2)).T
data2_2 = gaussian_filter(zvals_2, sigma_beta)

df_3= pd.DataFrame(dict(x3=dP50_vals_3, y3=f_stem_vals_3, z3=beta_3))
xcol_3, ycol_3, zcol_3 = "x3", "y3", "z3"
df_3 = df_3.sort_values(by=[xcol_3, ycol_3])
xvals_3 = df_3[xcol_3].unique()
yvals_3 = df_3[ycol_3].unique()
zvals_3 = df_3[zcol_3].values.reshape(len(xvals_3), len(yvals_3)).T
data2_3 = gaussian_filter(zvals_3, sigma_beta)

##Plot
fig, axs = plt.subplots(1, 3, figsize=(9,4), sharey=True)

CS1= axs[0].contourf(xvals_1, yvals_1, data2_1, levels_beta, cmap=cm_Ps, vmin=beta_min, vmax=beta_max)
CS1_1=axs[0].contour(xvals_1, yvals_1, data2_1, levels=beta_contours, colors='k')
axs[0].clabel(CS1_1, inline=1, fontsize=10, fmt='%.2f')
axs[0].set_ylim(min(yvals_1), max(yvals_1))
axs[0].set_title(r'Mean $P_{50}$ = -2.0 MPa')

CS2= axs[1].contourf(xvals_2, yvals_2, data2_2, levels_beta, cmap=cm_Ps, vmin=beta_min, vmax=beta_max)
CS2_1=axs[1].contour(xvals_2, yvals_2, data2_2, levels=beta_contours, colors='k')
axs[1].clabel(CS2_1, inline=1, fontsize=10, fmt='%.2f')
axs[1].set_ylim(min(yvals_2), max(yvals_2))
axs[1].set_title(r'Mean $P_{50}$ = -4.0 MPa')


CS3= axs[2].contourf(xvals_3, yvals_3, data2_3, levels_beta, cmap=cm_Ps, vmin=beta_min, vmax=beta_max)
CS3_1=axs[2].contour(xvals_3, yvals_3, data2_3, levels=beta_contours, colors='k')
axs[2].clabel(CS3_1, inline=1, fontsize=10, fmt='%.2f')
axs[2].set_ylim(min(yvals_3), max(yvals_3))
axs[2].set_title(r'Mean $P_{50}$ = -6.0 MPa')
CS3_map = ScalarMappable(norm=CS3.norm, cmap=CS3.cmap)

for ax in axs.flat:
    ax.set(xlabel=r'$\Delta P_{50}$ (MPa)', ylabel=r'$f_{stem,0}$')

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

plt.savefig('Other_Plots/experiment_2_beta_total.pdf')
plt.show()





