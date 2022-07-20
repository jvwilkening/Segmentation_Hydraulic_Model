from model_functions import import_traits_data, initialize_plant, ss_simulate_psi_k_norefilling, gen_k_from_P50
import numpy as np
from tqdm import tqdm
import scipy.stats as st
import os
from scipy.interpolate import interp1d

# set species traits, initialize plants with base traits from traits file
traits = import_traits_data()
plant = initialize_plant('GENERIC', traits)

# Define range of values for parameters that will be changed
num_mean_P50_vals = 30
num_delta_P50_vals = 20
sims_per_point = 50
mean_P50_min = -6.0
mean_P50_max = -2.0
dp50_spread = 1.5 #greatest delta P50 magnitude to use
f_stem = 0.75

VPD = 3.0


dp50_range = np.linspace(-dp50_spread, dp50_spread, num_delta_P50_vals)
mean_P50_vals = np.linspace(mean_P50_min, mean_P50_max, num_mean_P50_vals)

directory = "varying_segmentation_n_{0}_fstem_{1}".format(int(sims_per_point), int(f_stem*100.0))

#total_sets =num_stem_P50_vals * num_leaf_P50_vals
total_sets =num_mean_P50_vals * num_delta_P50_vals
gmin_frac = [.05]

k_actual_trial = np.zeros((sims_per_point, total_sets))
k_stem_actual_trial = np.zeros((sims_per_point, total_sets))
k_actual_frac_trial = np.zeros((sims_per_point, total_sets))
k_stem_frac_trial = np.zeros((sims_per_point, total_sets))
Px_trial = np.zeros((sims_per_point, total_sets))
P_soil_50_stem_trial = np.zeros((sims_per_point, total_sets))
P_soil_50_total_trial = np.zeros((sims_per_point, total_sets))
stem_P50_list=np.zeros(total_sets)
leaf_P50_list=np.zeros(total_sets)
mean_P50_list=np.zeros(total_sets)
delta_P50_list=np.zeros(total_sets)

set_number=0

for p50_set_num in tqdm(list(range(num_mean_P50_vals)), desc="p50_sets"):
    #plant.stem.P50_stem = stem_P50_vals[p50_set_num]
    #Pstem = [plant.stem.P50_stem]
    mean_P50 = mean_P50_vals[p50_set_num]
    g_frac = gmin_frac[0]
    gmax = plant.canopy.Gs_leaf
    #k_max_vals, a_vals = gen_k_from_P50(Pstem, plant.canopy.A_canopy, plant.canopy.H_canopy,
                                             #   sims_per_point)
    plant.stem.f_stem = f_stem
    plant.leaf.f_leaf = 1.0 - f_stem
    for delta_p50_set_num in tqdm(list(range(num_delta_P50_vals))):
        plant.stem.P50_stem = mean_P50 - dp50_range[delta_p50_set_num]/2.0
        Pstem = [plant.stem.P50_stem]
        plant.leaf.P50_leaf = mean_P50 + dp50_range[delta_p50_set_num]/2.0
        k_max_vals, a_vals = gen_k_from_P50(Pstem, plant.canopy.A_canopy, plant.canopy.H_canopy,
           sims_per_point)
        stem_P50_list[set_number] = plant.stem.P50_stem
        leaf_P50_list[set_number] = plant.leaf.P50_leaf
        mean_P50_list[set_number] = mean_P50
        delta_P50_list[set_number] = dp50_range[delta_p50_set_num]
        start_Psoil = 0.9 * max(Pstem[0], plant.leaf.P50_leaf)
        end_Psoil = min(Pstem[0], plant.leaf.P50_leaf)
        psi_soil_sets = max(int(abs(dp50_range[delta_p50_set_num]) * 10), 3)
        Psoil_range1 = np.linspace(start_Psoil, end_Psoil, psi_soil_sets)
        for k_sim_num in range(sims_per_point):
            k_max = k_max_vals[k_sim_num]
            a = a_vals[k_sim_num]
            plant.stem.k_stem_max = k_max / f_stem
            plant.leaf.k_leaf_max = k_max / (1.0 - f_stem)
            plant.stem.a_stem = a
            plant.leaf.a_leaf = a
            plant.canopy.Gs_min = g_frac * plant.canopy.Gs_leaf

            max_P50 = max(plant.stem.P50_stem, plant.leaf.P50_leaf)


            ps = (plant.stem.P50_stem + plant.leaf.P50_leaf)/2.0

            k_stem, k_actual, k_stem_frac, k_actual_frac, k_tissue_change, Px = ss_simulate_psi_k_norefilling(
                ps, plant, VPD)

            j = k_sim_num
            k = set_number

            k_actual_trial[j, k] = k_actual
            k_stem_actual_trial[j, k] = k_stem
            k_actual_frac_trial[j, k] = k_actual_frac
            k_stem_frac_trial[j, k] = k_stem_frac
            Px_trial[j, k] = Px

            stem_frac_trials = np.zeros(psi_soil_sets)
            total_frac_trials = np.zeros(psi_soil_sets)
            for psi_set in range(psi_soil_sets):
                psi_soil = Psoil_range1[psi_set]
                k_stem, k_actual, k_stem_frac, k_actual_frac, _, Px = ss_simulate_psi_k_norefilling(
                    psi_soil, plant, VPD)
                stem_frac_trials[psi_set] = k_stem_frac
                total_frac_trials[psi_set] = k_actual_frac

            if max(total_frac_trials) < 0.25:
                stem_frac_trials2 = np.zeros(psi_soil_sets)
                total_frac_trials2 = np.zeros(psi_soil_sets)
                start_Psoil = 0.5 * max(Pstem[0], plant.leaf.P50_leaf)
                end_Psoil = 0.85 * max(Pstem[0], plant.leaf.P50_leaf)
                Psoil_range2 = np.linspace(start_Psoil, end_Psoil, psi_soil_sets)
                for psi_set in range(psi_soil_sets):
                    psi_soil = Psoil_range2[psi_set]
                    k_stem, k_actual, k_stem_frac, k_actual_frac, _, Px = ss_simulate_psi_k_norefilling(
                        psi_soil, plant, VPD)
                    stem_frac_trials2[psi_set] = k_stem_frac
                    total_frac_trials2[psi_set] = k_actual_frac
                stem_frac_trials = np.append(stem_frac_trials, stem_frac_trials2)
                total_frac_trials = np.append(total_frac_trials, total_frac_trials2)
                Psoil_range = np.append(Psoil_range1, Psoil_range2)
            else:
                Psoil_range = Psoil_range1

            stem_func = interp1d(stem_frac_trials, Psoil_range, fill_value="extrapolate")
            total_func = interp1d(total_frac_trials, Psoil_range, fill_value="extrapolate")


            P_soil_50_stem_trial[j, k] = stem_func(0.5)
            P_soil_50_total_trial[j, k] = total_func(0.5)
        set_number = set_number + 1

k_nans = np.isnan(k_actual_trial)
k_actual_trial[k_nans] = 0.0
k_frac_nans = np.isnan(k_actual_frac_trial)
k_actual_frac_trial[k_frac_nans] = 0.0

mean_k = np.mean(k_actual_trial, axis=0)
std_error = st.sem(k_actual_trial, axis=0)
h = std_error * st.t.ppf((1.0+.95)/2.0, sims_per_point-1.0)

mean_stem_k = np.mean(k_stem_actual_trial, axis=0)
stem_std_error = st.sem(k_stem_actual_trial, axis=0)
h_stem = stem_std_error * st.t.ppf((1.0 + .95) / 2.0, sims_per_point - 1.0)


mean_k_frac = np.mean(k_actual_frac_trial, axis=0)
std_error_k_frac = st.sem(k_actual_frac_trial, axis=0)
h_k_frac = std_error_k_frac * st.t.ppf((1.0 + .95) / 2.0, sims_per_point - 1.0)
mean_stem_k_frac = np.mean(k_stem_frac_trial, axis=0)
stem_std_error_k_stem_frac = st.sem(k_stem_frac_trial, axis=0)
h_stem_frac = stem_std_error_k_stem_frac * st.t.ppf((1.0 + .95) / 2.0, sims_per_point - 1.0)

mean_stem50_Psoil = np.mean(P_soil_50_stem_trial, axis=0)
std_error_stem50_Psoil = st.sem(P_soil_50_stem_trial, axis=0)
h_stem50_Ps = std_error_stem50_Psoil * st.t.ppf((1.0 + .95) / 2.0, sims_per_point - 1.0)
mean_total50_Psoil = np.mean(P_soil_50_total_trial, axis=0)
std_error_total50_Psoil = st.sem(P_soil_50_total_trial, axis=0)
h_total0_Ps = std_error_total50_Psoil * st.t.ppf((1.0 + .95) / 2.0, sims_per_point - 1.0)

parent_dir = "Output"
path = os.path.join(parent_dir, directory)
os.makedirs(path)

np.savetxt(os.path.join(path, "k_plant_mean.txt"), mean_k, delimiter=',')
np.savetxt(os.path.join(path, "k_stem_mean.txt"), mean_stem_k, delimiter=',')
np.savetxt(os.path.join(path, "k_plant_ci.txt"), h, delimiter=',')
np.savetxt(os.path.join(path, "k_stem_ci.txt"), h_stem, delimiter=',')
np.savetxt(os.path.join(path, "k_plant_frac_mean.txt"), mean_k_frac, delimiter=',')
np.savetxt(os.path.join(path, "k_stem_frac_mean.txt"), mean_stem_k_frac, delimiter=',')
np.savetxt(os.path.join(path, "k_plant_frac_ci.txt"), h_k_frac, delimiter=',')
np.savetxt(os.path.join(path, "k_stem_frac_ci.txt"), h_stem_frac, delimiter=',')
np.savetxt(os.path.join(path, "stem_P50_vals.txt"), stem_P50_list, delimiter=',')
np.savetxt(os.path.join(path, "leaf_P50_vals.txt"), leaf_P50_list, delimiter=',')
np.savetxt(os.path.join(path, "stem50_Ps.txt"), mean_stem50_Psoil, delimiter=',')
np.savetxt(os.path.join(path, "total50_Ps.txt"), mean_total50_Psoil, delimiter=',')


