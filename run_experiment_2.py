from model_functions import import_traits_data, initialize_plant, ss_simulate_psi_k_norefilling, gen_k_from_P50
import numpy as np
from tqdm import tqdm
import scipy.stats as st
import os

#Where to save output
parent_dir = "Output"

# Define range of values for parameters that will be changed
fstem_param_vals= 50 #number of fstem values to run
safety_param_vals=50 #number of mean P50 vals to run
sims_per_point = 50 #simulations for each combination of fstem and P50
mean_P50_min = -6.0 #minimum mean P50 val [MPa]
mean_P50_max = -2.0
f_stem_min = 0.1
f_stem_max = 0.9
magnitude_dp50 = 0.1 #delta P50 magnitude [MPa]

#VPD in kPa
VPD=3.0


###Should not need to change anything below here

# set species traits, initialize plants with base traits from traits file
traits = import_traits_data()
conv_tree = initialize_plant('GENERIC', traits)
rev_tree = initialize_plant('GENERIC', traits)


f_stem_range = np.linspace(f_stem_min, f_stem_max, num=fstem_param_vals)
most_vulnerable_max = (2.0*mean_P50_max + magnitude_dp50)/2.0
most_vulnerable_min = (2.0 * mean_P50_min + magnitude_dp50) / 2.0
most_vulnerable_range = np.linspace(most_vulnerable_min, most_vulnerable_max, num=safety_param_vals)



def plot_ss_psi_s_v_k_MCs(plant, plant2, VPD, f_stem_range, most_vulnerable_range, fstem_param_vals, safety_param_vals, sims_per_point, magnitude_dp50, append=False):
    dp50 = -magnitude_dp50  # uses delta P50 as negative
    gmin_frac = [.05] #gsmin to gsmax ratio
    directory = "P50_vs_fstem_n_{0}_dp50_{1}".format(int(sims_per_point), int(-dp50*100.0))

    total_sets = fstem_param_vals * safety_param_vals

    k_conv_actual_trial = np.zeros((sims_per_point, total_sets))
    k_rev_actual_trial = np.zeros((sims_per_point, total_sets))
    k_conv_stem_actual_trial = np.zeros((sims_per_point, total_sets))
    k_rev_stem_actual_trial = np.zeros((sims_per_point, total_sets))
    k_conv_actual_frac_trial = np.zeros((sims_per_point, total_sets))
    k_rev_actual_frac_trial = np.zeros((sims_per_point, total_sets))
    k_conv_stem_frac_trial = np.zeros((sims_per_point, total_sets))
    k_rev_stem_frac_trial = np.zeros((sims_per_point, total_sets))
    Px_conv_trial = np.zeros((sims_per_point, total_sets))
    Px_rev_trial = np.zeros((sims_per_point, total_sets))
    max_P50_val=np.zeros(total_sets)
    fstem_val=np.zeros(total_sets)


    set_number=0

    for p50_set_num in tqdm(list(range(safety_param_vals)), desc="p50_sets"):
        max_P50 = most_vulnerable_range[p50_set_num] #most vulnerable P50 in the plant
        #max_P50 = mean_P50_range[p50_set_num]   # most vulnerable P50 in the plant
        g_frac = gmin_frac[0]
        gmax = plant.canopy.Gs_leaf
        conv_Pstem = [max_P50 + dp50]
        rev_Pstem = [max_P50]
        k_conv_max_vals, a_conv_vals = gen_k_from_P50(conv_Pstem, plant.canopy.A_canopy, plant.canopy.H_canopy,
                                                    sims_per_point)
        k_rev_max_vals, a_rev_vals = gen_k_from_P50(rev_Pstem, plant.canopy.A_canopy, plant.canopy.H_canopy,
                                                    sims_per_point)
        for f_stem_set_num in tqdm(list(range(fstem_param_vals)), desc="f_sets"):
            f_stem = f_stem_range[f_stem_set_num]
            max_P50_val[set_number] = max_P50
            fstem_val[set_number] = f_stem
            for k_sim_num in range(sims_per_point):

                k_conv_max = k_conv_max_vals[k_sim_num]
                a_conv = a_conv_vals[k_sim_num]
                k_rev_max = k_rev_max_vals[k_sim_num]
                a_rev = a_rev_vals[k_sim_num]

                plant.stem.f_stem = f_stem
                plant.leaf.f_leaf = 1.0 - f_stem
                plant.stem.k_stem_max = k_conv_max / f_stem
                plant.leaf.k_leaf_max = k_conv_max / (1.0 - f_stem)
                plant.stem.a_stem = a_conv
                plant.leaf.a_leaf = a_conv
                plant.stem.P50_stem = conv_Pstem[0]
                plant.leaf.P50_leaf = conv_Pstem[0] - dp50
                plant.canopy.Gs_min = g_frac * plant.canopy.Gs_leaf

                plant2.stem.f_stem = f_stem
                plant2.leaf.f_leaf = 1.0 - f_stem
                plant2.stem.k_stem_max = k_rev_max / f_stem
                plant2.leaf.k_leaf_max = k_rev_max / (1.0 - f_stem)
                plant2.stem.a_stem = a_rev
                plant2.leaf.a_leaf = a_rev
                plant2.stem.P50_stem = rev_Pstem[0]
                plant2.leaf.P50_leaf = rev_Pstem[0] + dp50
                plant2.canopy.Gs_min = g_frac * plant2.canopy.Gs_leaf

                k_conv_stem, k_conv_actual, k_conv_stem_frac, k_conv_actual_frac, _, Pxconv = ss_simulate_psi_k_norefilling(
                    max_P50, plant, VPD)
                k_rev_stem, k_rev_actual, k_rev_stem_frac, k_rev_actual_frac, _, Pxrev = ss_simulate_psi_k_norefilling(
                    max_P50, plant2, VPD)

                j = k_sim_num
                k = set_number

                k_conv_actual_trial[j, k] = k_conv_actual
                k_conv_stem_actual_trial[j, k] = k_conv_stem
                k_conv_actual_frac_trial[j, k] = k_conv_actual_frac
                k_conv_stem_frac_trial[j, k] = k_conv_stem_frac
                k_rev_actual_trial[j, k] = k_rev_actual
                k_rev_stem_actual_trial[j, k] = k_rev_stem
                k_rev_actual_frac_trial[j, k] = k_rev_actual_frac
                k_rev_stem_frac_trial[j, k] = k_rev_stem_frac
                Px_conv_trial[j, k] = Pxconv
                Px_rev_trial[j, k] = Pxrev
            set_number=set_number+1



    #replace nans that sometime rev up in dead plants in clay soils with 0.0
    k_conv_nans = np.isnan(k_conv_actual_trial)
    k_conv_actual_trial[k_conv_nans] = 0.0
    k_rev_nans = np.isnan(k_rev_actual_trial)
    k_rev_actual_trial[k_rev_nans] = 0.0
    k_conv_frac_nans = np.isnan(k_conv_actual_frac_trial)
    k_conv_actual_frac_trial[k_conv_frac_nans] = 0.0
    k_rev_frac_nans = np.isnan(k_rev_actual_frac_trial)
    k_rev_actual_frac_trial[k_rev_frac_nans] = 0.0

    #calculate statistics across all sims
    mean_conv_k = np.mean(k_conv_actual_trial, axis=0)
    mean_rev_k = np.mean(k_rev_actual_trial, axis=0)
    conv_std_error = st.sem(k_conv_actual_trial, axis=0)
    rev_std_error = st.sem(k_rev_actual_trial, axis=0)
    h_conv = conv_std_error * st.t.ppf((1.0+.95)/2.0, sims_per_point-1.0)
    h_rev = rev_std_error * st.t.ppf((1.0 + .95) / 2.0, sims_per_point - 1.0)
    mean_conv_stem_k = np.mean(k_conv_stem_actual_trial, axis=0)
    mean_rev_stem_k = np.mean(k_rev_stem_actual_trial, axis=0)
    conv_stem_std_error = st.sem(k_conv_stem_actual_trial, axis=0)
    rev_stem_std_error = st.sem(k_rev_stem_actual_trial, axis=0)
    h_conv_stem = conv_stem_std_error * st.t.ppf((1.0 + .95) / 2.0, sims_per_point - 1.0)
    h_rev_stem = rev_stem_std_error * st.t.ppf((1.0 + .95) / 2.0, sims_per_point - 1.0)


    mean_conv_k_frac = np.mean(k_conv_actual_frac_trial, axis=0)
    mean_rev_k_frac = np.mean(k_rev_actual_frac_trial, axis=0)
    conv_std_error_k_frac = st.sem(k_conv_actual_frac_trial, axis=0)
    rev_std_error_k_frac = st.sem(k_rev_actual_frac_trial, axis=0)
    h_conv_k_frac = conv_std_error_k_frac * st.t.ppf((1.0 + .95) / 2.0, sims_per_point - 1.0)
    h_rev_k_frac = rev_std_error_k_frac * st.t.ppf((1.0 + .95) / 2.0, sims_per_point - 1.0)
    mean_conv_stem_k_frac = np.mean(k_conv_stem_frac_trial, axis=0)
    mean_rev_stem_k_frac = np.mean(k_rev_stem_frac_trial, axis=0)
    conv_stem_std_error_k_stem_frac = st.sem(k_conv_stem_frac_trial, axis=0)
    rev_stem_std_error_k_stem_frac = st.sem(k_rev_stem_frac_trial, axis=0)
    h_conv_stem_frac = conv_stem_std_error_k_stem_frac * st.t.ppf((1.0 + .95) / 2.0, sims_per_point - 1.0)
    h_rev_stem_frac = rev_stem_std_error_k_stem_frac * st.t.ppf((1.0 + .95) / 2.0, sims_per_point - 1.0)


    #Save output files

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
    np.savetxt(os.path.join(path, "rev_k_plant_frac_mean.txt"), mean_rev_k_frac, delimiter=',')
    np.savetxt(os.path.join(path, "rev_k_stem_frac_mean.txt"), mean_rev_stem_k_frac, delimiter=',')
    np.savetxt(os.path.join(path, "rev_k_plant_frac_ci.txt"), h_rev_k_frac, delimiter=',')
    np.savetxt(os.path.join(path, "rev_k_stem_frac_ci.txt"), h_rev_stem_frac, delimiter=',')
    np.savetxt(os.path.join(path, "max_P50_vals.txt"), max_P50_val, delimiter=',')
    np.savetxt(os.path.join(path, "f_stem_vals.txt"), fstem_val, delimiter=',')


if __name__ == '__main__':
    plot_ss_psi_s_v_k_MCs(conv_tree, rev_tree, VPD, f_stem_range, most_vulnerable_range, fstem_param_vals, safety_param_vals, sims_per_point, magnitude_dp50)

