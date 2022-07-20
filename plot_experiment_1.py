#plots results in terms of lambda stem and lamda soil vs psi soil for Experiment I

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import scipy.stats as st

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42


#Define arameter values used in simulations to plot (used to retrieve results files)

f_vals = [25, 50, 75] #fstem vals used in simulations (input as % not fraction)
dp50 = 1.0 # delta P50 val used
p_vals = [2.5, 4.5, 6.5] #max P50 vals used in simulations
n = 750.0 #number of sims

plot_cols=[plt.cm.cool(c) for c in np.linspace(0,1,len(f_vals))]
fig, axs = plt.subplots(1, 3, sharey=True, sharex=True, figsize=(8,5))

##plot lambda total
for k in range(len(f_vals)):
    f=f_vals[k]
    for j in range(len(p_vals)):
        p=p_vals[j]
        filename = "Output/n_{0}_p50_{1}_fstem_{2}_dp50{3}".format(int(n),int(p), int(f), int(dp50))
        rev_stem = np.loadtxt('%s/rev_k_stem_frac_mean.txt' % (filename))
        conv_stem = np.loadtxt('%s/conv_k_stem_frac_mean.txt' % (filename))
        conv_stem_ci = np.loadtxt('%s/conv_k_stem_frac_ci.txt' % (filename))
        rev_stem_ci = np.loadtxt('%s/rev_k_stem_frac_ci.txt' % (filename))
        rev_plant = np.loadtxt('%s/rev_k_plant_frac_mean.txt' % (filename))
        conv_plant = np.loadtxt('%s/conv_k_plant_frac_mean.txt' % (filename))
        conv_plant_ci = np.loadtxt('%s/conv_k_plant_frac_ci.txt' % (filename))
        rev_plant_ci = np.loadtxt('%s/rev_k_plant_frac_ci.txt' % (filename))
        s = np.loadtxt('%s/s_range.txt' % (filename))
        num_samples = n


        plant_diff = np.zeros_like(rev_stem)
        mean_diff_CI = np.zeros_like(rev_stem)
        s_adj = np.zeros_like(rev_stem)
        for l in range(len(rev_plant)):
            sd_rev = rev_plant_ci[l]/st.t.ppf((1.0 + .95) / 2.0, num_samples - 1.0)*np.sqrt(num_samples)
            sd_conv = conv_plant_ci[l]/st.t.ppf((1.0 + .95) / 2.0, num_samples - 1.0)*np.sqrt(num_samples)
            pooled_sd = np.sqrt((((num_samples-1.0)*sd_rev**2.0) + ((num_samples-1.0)*sd_conv**2.0))/(num_samples+num_samples-2.0))
            mean_diff_CI[l] = 1.96*pooled_sd*np.sqrt((1.0/num_samples)+(1.0/num_samples))
            plant_diff[l] = conv_plant[l] - rev_plant[l]
            s_adj[l] = s[l] + p
        axs[k].plot(s_adj, plant_diff, color=plot_cols[j], label= 'f_stem = %0.2f' % (f_vals[k]/100.0))
        axs[k].fill_between(s_adj, plant_diff+mean_diff_CI, plant_diff-mean_diff_CI, color=plot_cols[j], alpha=0.2)

for ax in axs.flat:
    ax.set(xlabel=r'$\psi_{soil}$ relative to Mean $P_{50}$ [MPa]', ylabel=r'$\lambda_{total}$ (%)')
    ax.axhline(y=0.0, color='gray', alpha=0.5)
    ax.set_xticks([-2.0, -1.0, 0.0, 1.0, 2.0])
    ax.set_xlim(-1.5, 1.5)
for ax in axs.flat: #only show outer axis labels
    ax.label_outer()

legend_elements = [Line2D([0], [0], color=plot_cols[0], label= 'P50 = -%2.1f MPa' % (p_vals[0]),
                          lw=3),
                   Line2D([0], [0], color=plot_cols[1], label= 'P50 = -%2.1f MPa' % (p_vals[1]),
                          lw=3),
                   Line2D([0], [0], color=plot_cols[2], label= 'P50 = -%2.1f MPa' % p_vals[2],
                          lw=3),
                   ]
fig.legend(handles=legend_elements, loc="upper left")

#Update labels as necessary
axs[0].set_title('f_stem,0=0.25')
axs[1].set_title('f_stem,0=0.5')
axs[2].set_title('f_stem,0=0.75')

pad = 5
fig.subplots_adjust(left=0.22, top=0.88)
plt.savefig('Results_Plots/experiment_1_lambda_total.pdf')
plt.show()


##plot lambda stem
plot_cols=[plt.cm.cool(c) for c in np.linspace(0,1,len(f_vals))]
fig, axs = plt.subplots(1, 3, sharey=True, sharex=True, figsize=(8,5))


for k in range(len(f_vals)):
    f=f_vals[k]
    for j in range(len(p_vals)):
        p=p_vals[j]
        filename = "Output/n_{0}_p50_{1}_fstem_{2}_dp50{3}".format(int(n),int(p), int(f), int(dp50))
        rev_stem = np.loadtxt('%s/rev_k_stem_frac_mean.txt' % (filename))*100.0
        conv_stem = np.loadtxt('%s/conv_k_stem_frac_mean.txt' % (filename))*100.0
        conv_stem_ci = np.loadtxt('%s/conv_k_stem_frac_ci.txt' % (filename))*100.0
        rev_stem_ci = np.loadtxt('%s/rev_k_stem_frac_ci.txt' % (filename))*100.0
        rev_plant = np.loadtxt('%s/rev_k_plant_frac_mean.txt' % (filename))
        conv_plant = np.loadtxt('%s/conv_k_plant_frac_mean.txt' % (filename))
        conv_plant_ci = np.loadtxt('%s/conv_k_plant_frac_ci.txt' % (filename))
        rev_plant_ci = np.loadtxt('%s/rev_k_plant_frac_ci.txt' % (filename))
        s = np.loadtxt('%s/s_range.txt' % (filename))
        num_samples = n


        stem_diff = np.zeros_like(rev_stem)
        mean_diff_CI = np.zeros_like(rev_stem)
        s_adj = np.zeros_like(rev_stem)
        for l in range(len(rev_stem)):
            sd_rev = rev_stem_ci[l]/st.t.ppf((1.0 + .95) / 2.0, num_samples - 1.0)*np.sqrt(num_samples)
            sd_conv = conv_stem_ci[l]/st.t.ppf((1.0 + .95) / 2.0, num_samples - 1.0)*np.sqrt(num_samples)
            pooled_sd = np.sqrt((((num_samples-1.0)*sd_rev**2.0) + ((num_samples-1.0)*sd_conv**2.0))/(num_samples+num_samples-2.0))
            mean_diff_CI[l] = 1.96*pooled_sd*np.sqrt((1.0/num_samples)+(1.0/num_samples))
            stem_diff[l] = conv_stem[l] - rev_stem[l]
            s_adj[l] = s[l] + p
        axs[k].plot(s_adj, stem_diff, color=plot_cols[j], label= 'f_stem = %0.2f' % (f_vals[k]/100.0))
        axs[k].fill_between(s_adj, stem_diff+mean_diff_CI, stem_diff-mean_diff_CI, color=plot_cols[j], alpha=0.2)

for ax in axs.flat:
    ax.set(xlabel=r'$\psi_{soil}$ relative to Mean $P_{50}$ [MPa]', ylabel=r'$\lambda_{stem}$ (%)')
    ax.axhline(y=0.0, color='gray', alpha=0.5)
    ax.set_xticks([-2.0, -1.0, 0.0, 1.0, 2.0])
    ax.set_xlim(-1.5, 1.5)
for ax in axs.flat: #only show outer axis labels
    ax.label_outer()

legend_elements = [Line2D([0], [0], color=plot_cols[0], label= 'P50 = -%2.1f MPa' % (p_vals[0]),
                          lw=3),
                   Line2D([0], [0], color=plot_cols[1], label= 'P50 = -%2.1f MPa' % (p_vals[1]),
                          lw=3),
                   Line2D([0], [0], color=plot_cols[2], label= 'P50 = -%2.1f MPa' % p_vals[2],
                          lw=3),
                   ]
fig.legend(handles=legend_elements, loc="upper left")

#update labels as necessary
axs[0].set_title('f_stem,0=0.25')
axs[1].set_title('f_stem,0=0.5')
axs[2].set_title('f_stem,0=0.75')

pad = 5 # in points


fig.subplots_adjust(left=0.22, top=0.88)
plt.savefig('Results_Plots/experiment_1_lambda_stem.pdf')
plt.show()


