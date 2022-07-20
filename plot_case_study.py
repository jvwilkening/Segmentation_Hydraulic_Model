import matplotlib.pyplot as plt
import numpy as np
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42


###Update all file names as necessary
oak_total_25_obs = np.loadtxt('oak_obs_fstem_25/k_plant_frac_mean.txt')
oak_total_50_obs = np.loadtxt('oak_obs_fstem_50/k_plant_frac_mean.txt')
oak_total_75_obs = np.loadtxt('oak_obs_fstem_75/k_plant_frac_mean.txt')
psi_s = np.loadtxt('oak_obs_fstem_25/psi_s.txt')

oak_stem_25_obs = np.loadtxt('oak_obs_fstem_25/k_stem_frac_mean.txt')
oak_stem_50_obs = np.loadtxt('oak_obs_fstem_50/k_stem_frac_mean.txt')
oak_stem_75_obs = np.loadtxt('oak_obs_fstem_75/k_stem_frac_mean.txt')

oak_total_25_rev = np.loadtxt('oak_rev_fstem_25/k_plant_frac_mean.txt')
oak_total_50_rev = np.loadtxt('oak_rev_fstem_50/k_plant_frac_mean.txt')
oak_total_75_rev = np.loadtxt('oak_rev_fstem_75/k_plant_frac_mean.txt')

oak_stem_25_rev = np.loadtxt('oak_rev_fstem_25/k_stem_frac_mean.txt')
oak_stem_50_rev = np.loadtxt('oak_rev_fstem_50/k_stem_frac_mean.txt')
oak_stem_75_rev = np.loadtxt('oak_rev_fstem_75/k_stem_frac_mean.txt')

pop_total_25_obs = np.loadtxt('pop_obs_fstem_25/k_plant_frac_mean.txt')
pop_total_50_obs = np.loadtxt('pop_obs_fstem_50/k_plant_frac_mean.txt')
pop_total_75_obs = np.loadtxt('pop_obs_fstem_75/k_plant_frac_mean.txt')
pop_psi_s = np.loadtxt('pop_obs_fstem_25/psi_s.txt')

pop_stem_25_obs = np.loadtxt('pop_obs_fstem_25/k_stem_frac_mean.txt')
pop_stem_50_obs = np.loadtxt('pop_obs_fstem_50/k_stem_frac_mean.txt')
pop_stem_75_obs = np.loadtxt('pop_obs_fstem_75/k_stem_frac_mean.txt')

pop_total_25_rev = np.loadtxt('pop_rev_fstem_25/k_plant_frac_mean.txt')
pop_total_50_rev = np.loadtxt('pop_rev_fstem_50/k_plant_frac_mean.txt')
pop_total_75_rev = np.loadtxt('pop_rev_fstem_75/k_plant_frac_mean.txt')

pop_stem_25_rev = np.loadtxt('pop_rev_fstem_25/k_stem_frac_mean.txt')
pop_stem_50_rev = np.loadtxt('pop_rev_fstem_50/k_stem_frac_mean.txt')
pop_stem_75_rev = np.loadtxt('pop_rev_fstem_75/k_stem_frac_mean.txt')

fig, axs = plt.subplots(3, 2, figsize=(8,9))


axs[2,0].plot(psi_s, pop_total_75_obs, lw=3, color='deepskyblue', alpha=1.0, label="Total Observed", linewidth=2)
axs[2,0].plot(psi_s, pop_total_75_rev, lw=3, color='aqua', alpha=1.0, label="Total Flipped", linestyle='dotted')
axs[2,0].plot(psi_s, pop_stem_75_obs, lw=3, color='mediumvioletred', alpha=1.0, label="Stem Observed", linewidth=2)
axs[2,0].plot(psi_s, pop_stem_75_rev, lw=3, color='hotpink', alpha=1.0, label="Stem Flipped", linestyle='dotted')
axs[2,0].set_xlim((-4.0, 0.0))


axs[0,0].plot(psi_s, pop_total_25_obs, lw=3, color='deepskyblue', alpha=1.0, label="Total Observed", linewidth=2)
axs[0,0].plot(psi_s, pop_total_25_rev, lw=3, color='aqua', alpha=1.0, label="Total Flipped", linestyle='dotted')
axs[0,0].plot(psi_s, pop_stem_25_obs, lw=3, color='mediumvioletred', alpha=1.0, label="Stem Observed", linewidth=2)
axs[0,0].plot(psi_s, pop_stem_25_rev, lw=3, color='hotpink', alpha=1.0, label="Stem Flipped", linestyle='dotted')
axs[0,0].set_xlim((-4.0, 0.0))


axs[1,0].plot(psi_s, pop_total_50_obs, lw=3, color='deepskyblue', alpha=1.0, label="Total Observed", linewidth=2)
axs[1,0].plot(psi_s, pop_total_50_rev, lw=3, color='aqua', alpha=1.0, label="Total Flipped", linestyle='dotted')
axs[1,0].plot(psi_s, pop_stem_50_obs, lw=3, color='mediumvioletred', alpha=1.0, label="Stem Observed", linewidth=2)
axs[1,0].plot(psi_s, pop_stem_50_rev, lw=3, color='hotpink', alpha=1.0, label="Stem Flipped", linestyle='dotted')
axs[1,0].set_xlim((-4.0, 0.0))



axs[2,1].plot(psi_s, oak_total_75_obs, lw=3, color='deepskyblue', alpha=1.0, label="Total Observed", linewidth=2)
axs[2,1].plot(psi_s, oak_total_75_rev, lw=3, color='aqua', alpha=1.0, label="Total Flipped", linestyle='dotted')
axs[2,1].plot(psi_s, oak_stem_75_obs, lw=3, color='mediumvioletred', alpha=1.0, label="Stem Observed", linewidth=2)
axs[2,1].plot(psi_s, oak_stem_75_rev, lw=3, color='hotpink', alpha=1.0, label="Stem Flipped", linestyle='dotted')
axs[2,1].set_xlim((-6.0, -2.0))


axs[0,1].plot(psi_s, oak_total_25_obs, lw=3, color='deepskyblue', alpha=1.0, label="Total Observed", linewidth=2)
axs[0,1].plot(psi_s, oak_total_25_rev, lw=3, color='aqua', alpha=1.0, label="Total Flipped", linestyle='dotted')
axs[0,1].plot(psi_s, oak_stem_25_obs, lw=3, color='mediumvioletred', alpha=1.0, label="Stem Observed", linewidth=2)
axs[0,1].plot(psi_s, oak_stem_25_rev, lw=3, color='hotpink', alpha=1.0, label="Stem Flipped", linestyle='dotted')
axs[0,1].set_xlim((-6.0, -2.0))


axs[1,1].plot(psi_s, oak_total_50_obs, lw=3, color='deepskyblue', alpha=1.0, label="Total Observed", linewidth=2)
axs[1,1].plot(psi_s, oak_total_50_rev, lw=3, color='aqua', alpha=1.0, label="Total Flipped", linestyle='dotted')
axs[1,1].plot(psi_s, oak_stem_50_obs, lw=3, color='mediumvioletred', alpha=1.0, label="Stem Observed", linewidth=2)
axs[1,1].plot(psi_s, oak_stem_50_rev, lw=3, color='hotpink', alpha=1.0, label="Stem Flipped", linestyle='dotted')
axs[1,1].set_xlim((-6.0, -2.0))
for ax in axs.flat:
    ax.set(xlabel='Soil Water Potential (-MPa)', ylabel='Fraction Conductance Remaining')

for ax in axs.flat:
    ax.label_outer()
plt.savefig('../Other_Plots/pop_oak_v_soil_potential.pdf')
plt.show()
