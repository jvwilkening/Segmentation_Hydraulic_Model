from model_functions import import_traits_data, initialize_plant, ss_simulate_psi_k_norefilling, gen_k_from_P50
import numpy as np
from tqdm import tqdm
import scipy.stats as st
import os

#parent directory for output
parent_dir = "Output"

# Define values for independent parameters
most_vulnerable = [-2.5] #most vulnerable P50 value in plants
f_stem = 0.25 #fraction of total hydraulic resistance in stems, must be >0 and <1
delta_P50 = 0.0 #difference between P50s in tissues (conventional - reverse)

VPD=3.0 #VPD in kPa

# Define range of psi soil values to run simulations over
psi_s_min = most_vulnerable[0] - 2.0 #psi_soil range defaults to +/- 2 MPa around max P50
psi_s_max = np.minimum(0.0, (most_vulnerable[0] + 2.0))
num_s = 40 # number of points over s range
num_k_sims = 300 # conductance values to generate for each P50 value
s_range = np.linspace(psi_s_min, psi_s_max, num=num_s)

###Should not need to change anything below here

# set species traits, initialize plants with base traits from traits file [P50 difference set later]
traits = import_traits_data()
conv_tree = initialize_plant('GENERIC', traits)
rev_tree = initialize_plant('GENERIC', traits)
def plot_ss_psi_s_v_k_MCs(plant, plant2, VPD, s_range, num_k_sims):
    #directory where all output will be saved
    directory = "n_{0}_p50_{1}_fstem_{2}_dP50_{3}".format(int(num_k_sims), int(-most_vulnerable[0]), int(f_stem*100.0), int(delta_P50*10.0))

    dp50 = -delta_P50
    gmin_frac = [.05] #ratio of gsmin to gsmax
    conv_Pstem_range = [x + dp50 for x in most_vulnerable] #plant with more vulnerable leaves
    rev_Pstem_range = most_vulnerable #plant with more vulnerable stem
    #generate other plant parameters based on values of prescribed parameters
    k_conv_max_vals, a_conv_vals = gen_k_from_P50(conv_Pstem_range, plant.canopy.A_canopy, plant.canopy.H_canopy, num_k_sims)
    k_rev_max_vals, a_rev_vals = gen_k_from_P50(rev_Pstem_range, plant.canopy.A_canopy, plant.canopy.H_canopy, num_k_sims)


    for i in tqdm(list(range(len(most_vulnerable))), desc="f_sets"):
        g_frac = gmin_frac[i]
        k_conv_actual_trial = np.zeros((num_k_sims, num_s))
        k_rev_actual_trial = np.zeros((num_k_sims, num_s))
        k_conv_stem_actual_trial = np.zeros((num_k_sims, num_s))
        k_rev_stem_actual_trial = np.zeros((num_k_sims, num_s))
        k_conv_actual_frac_trial = np.zeros((num_k_sims, num_s))
        k_rev_actual_frac_trial = np.zeros((num_k_sims, num_s))
        k_conv_stem_frac_trial = np.zeros((num_k_sims, num_s))
        k_rev_stem_frac_trial = np.zeros((num_k_sims, num_s))
        Px_conv_trial = np.zeros((num_k_sims, num_s))
        Px_rev_trial = np.zeros((num_k_sims, num_s))
        rev_P50_stem=rev_Pstem_range[i]
        conv_P50_stem=conv_Pstem_range[i]
        for j in tqdm(list(range(num_k_sims)), desc="k_sims"):
            set_num = j + (i*num_k_sims)
            k_conv_max = k_conv_max_vals[set_num]
            a_conv = a_conv_vals[set_num]
            k_rev_max = k_rev_max_vals[set_num]
            a_rev = a_rev_vals[set_num]
            plant.stem.f_stem = f_stem
            plant.leaf.f_leaf = 1.0 - f_stem
            plant.stem.k_stem_max = k_conv_max / f_stem
            plant.leaf.k_leaf_max = k_conv_max / (1.0 - f_stem)
            plant.stem.a_stem = a_conv
            plant.leaf.a_leaf = a_conv
            plant.stem.P50_stem = conv_P50_stem
            plant.leaf.P50_leaf = conv_P50_stem - dp50
            plant.canopy.Gs_min = g_frac*plant.canopy.Gs_leaf

            plant2.stem.f_stem = f_stem
            plant2.leaf.f_leaf = 1.0 - f_stem
            plant2.stem.k_stem_max = k_rev_max / f_stem
            plant2.leaf.k_leaf_max = k_rev_max / (1.0 - f_stem)
            plant2.stem.a_stem = a_rev
            plant2.leaf.a_leaf = a_rev
            plant2.stem.P50_stem = rev_P50_stem
            plant2.leaf.P50_leaf = rev_P50_stem + dp50
            plant2.canopy.Gs_min = g_frac * plant2.canopy.Gs_leaf
            for k in range(len(s_range)):
                psi_s0 = s_range[k]
                k_conv_stem, k_conv_actual, k_conv_stem_frac, k_conv_actual_frac, _, Pxconv = ss_simulate_psi_k_norefilling(psi_s0, plant, VPD)
                k_rev_stem, k_rev_actual, k_rev_stem_frac, k_rev_actual_frac, _, _change_rev, Pxrev = ss_simulate_psi_k_norefilling(psi_s0, plant2, VPD)
                k_conv_actual_trial[j, k] = k_conv_actual
                k_conv_stem_actual_trial[j, k] = k_conv_stem
                k_conv_actual_frac_trial[j, k] = k_conv_actual_frac
                k_conv_stem_frac_trial[j, k] = k_conv_stem_frac
                k_rev_actual_trial[j, k] = k_rev_actual
                k_rev_stem_actual_trial[j, k] = k_rev_stem
                k_rev_actual_frac_trial[j, k] = k_rev_actual_frac
                k_rev_stem_frac_trial[j, k] = k_rev_stem_frac
                Px_conv_trial[j,k] = Pxconv
                Px_rev_trial[j,k] = Pxrev

        #replace nans that sometime rev up in dead plants in clay soils with 0.0
        k_conv_nans = np.isnan(k_conv_actual_trial)
        k_conv_actual_trial[k_conv_nans] = 0.0
        k_rev_nans = np.isnan(k_rev_actual_trial)
        k_rev_actual_trial[k_rev_nans] = 0.0
        k_conv_frac_nans = np.isnan(k_conv_actual_frac_trial)
        k_conv_actual_frac_trial[k_conv_frac_nans] = 0.0
        k_rev_frac_nans = np.isnan(k_rev_actual_frac_trial)
        k_rev_actual_frac_trial[k_rev_frac_nans] = 0.0

        #calculate statistics over all runs
        mean_conv_k = np.mean(k_conv_actual_trial, axis=0)
        mean_rev_k = np.mean(k_rev_actual_trial, axis=0)
        conv_std_error = st.sem(k_conv_actual_trial, axis=0)
        rev_std_error = st.sem(k_rev_actual_trial, axis=0)
        h_conv = conv_std_error * st.t.ppf((1.0+.95)/2.0, num_k_sims-1.0)
        h_rev = rev_std_error * st.t.ppf((1.0 + .95) / 2.0, num_k_sims - 1.0)
        mean_conv_stem_k = np.mean(k_conv_stem_actual_trial, axis=0)
        mean_rev_stem_k = np.mean(k_rev_stem_actual_trial, axis=0)
        conv_stem_std_error = st.sem(k_conv_stem_actual_trial, axis=0)
        rev_stem_std_error = st.sem(k_rev_stem_actual_trial, axis=0)
        h_conv_stem = conv_stem_std_error * st.t.ppf((1.0 + .95) / 2.0, num_k_sims - 1.0)
        h_rev_stem = rev_stem_std_error * st.t.ppf((1.0 + .95) / 2.0, num_k_sims - 1.0)

        mean_conv_k_frac = np.mean(k_conv_actual_frac_trial, axis=0)
        mean_rev_k_frac = np.mean(k_rev_actual_frac_trial, axis=0)
        conv_std_error_k_frac = st.sem(k_conv_actual_frac_trial, axis=0)
        rev_std_error_k_frac = st.sem(k_rev_actual_frac_trial, axis=0)
        h_conv_k_frac = conv_std_error_k_frac * st.t.ppf((1.0 + .95) / 2.0, num_k_sims - 1.0)
        h_rev_k_frac = rev_std_error_k_frac * st.t.ppf((1.0 + .95) / 2.0, num_k_sims - 1.0)
        mean_conv_stem_k_frac = np.mean(k_conv_stem_frac_trial, axis=0)
        mean_rev_stem_k_frac = np.mean(k_rev_stem_frac_trial, axis=0)
        conv_stem_std_error_k_stem_frac = st.sem(k_conv_stem_frac_trial, axis=0)
        rev_stem_std_error_k_stem_frac = st.sem(k_rev_stem_frac_trial, axis=0)
        h_conv_stem_frac = conv_stem_std_error_k_stem_frac * st.t.ppf((1.0 + .95) / 2.0, num_k_sims - 1.0)
        h_rev_stem_frac = rev_stem_std_error_k_stem_frac * st.t.ppf((1.0 + .95) / 2.0, num_k_sims - 1.0)

    path = os.path.join(parent_dir, directory)
    os.makedirs(path)

    np.savetxt(os.path.join(path, "conv_k_plant_mean.txt"), mean_conv_k, delimiter=',')
    np.savetxt(os.path.join(path, "conv_k_stem_mean.txt"), mean_conv_stem_k, delimiter=',')
    np.savetxt(os.path.join(path, "conv_k_plant_ci.txt"), h_conv, delimiter=',')
    np.savetxt(os.path.join(path, "conv_k_stem_ci.txt"), h_conv_stem, delimiter=',')
    np.savetxt(os.path.join(path, "conv_k_plant_frac_mean.txt"), mean_conv_k_frac, delimiter=',')
    np.savetxt(os.path.join(path, "conv_k_stem_frac_mean.txt"), mean_conv_stem_k_frac, delimiter=',')
    np.savetxt(os.path.join(path, "conv_k_plant_frac_ci.txt"), h_conv_k_frac, delimiter=',')
    np.savetxt(os.path.join(path, "conv_k_stem_frac_ci.txt"), h_conv_stem_frac, delimiter=',')


    np.savetxt(os.path.join(path, "rev_k_plant_mean.txt"), mean_rev_k, delimiter=',')
    np.savetxt(os.path.join(path, "rev_k_stem_mean.txt"), mean_rev_stem_k, delimiter=',')
    np.savetxt(os.path.join(path, "rev_k_plant_ci.txt"), h_rev, delimiter=',')
    np.savetxt(os.path.join(path, "rev_k_stem_ci.txt"), h_rev_stem, delimiter=',')
    np.savetxt(os.path.join(path, "s_range.txt"), s_range, delimiter=',')
    np.savetxt(os.path.join(path, "rev_k_plant_frac_mean.txt"), mean_rev_k_frac, delimiter=',')
    np.savetxt(os.path.join(path, "rev_k_stem_frac_mean.txt"), mean_rev_stem_k_frac, delimiter=',')
    np.savetxt(os.path.join(path, "rev_k_plant_frac_ci.txt"), h_rev_k_frac, delimiter=',')
    np.savetxt(os.path.join(path, "rev_k_stem_frac_ci.txt"), h_rev_stem_frac, delimiter=',')

if __name__ == '__main__':
    plot_ss_psi_s_v_k_MCs(conv_tree, rev_tree, VPD, s_range, num_k_sims)
