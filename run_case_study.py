# generates k vs soil curves for specific plants (pop and oak flipped vs observed)

from model_functions import import_traits_data, initialize_plant, ss_simulate_psi_k_norefilling
import numpy as np
from tqdm import tqdm
import os

# set species traits, initialize plants with base traits from traits file
#Update which tree and configuration (observed vs reverse) to run in function call at end
traits = import_traits_data()
oak = initialize_plant('OAK', traits)
pop = initialize_plant('POP', traits)
oak_rev = initialize_plant('OAK_REV', traits)
pop_rev = initialize_plant('POP_REV', traits)

# Define range of values for parameters that will be changed
f_range = [.25]  # fraction of resistance in stem, just do one at a time
VPD = 3.0 #VPD in kPa

# Define range of psi_soil and number of points
psi_s_min = -6.0
psi_s_max = 0.0
num_s = 50  # number of points over s range

s_range = np.linspace(psi_s_min, psi_s_max, num=num_s)

#update directory name as necessary
directory = "pop_obs_fstem_{0}".format(int(f_range[0] * 100.0))
parent_dir = "Output"
path = os.path.join(parent_dir, directory)
os.makedirs(path)


def plot_ss_soil_vs_k(base_plant, VPD, s_range):

    k_stem_max = base_plant.stem.k_stem_max
    k_leaf_max = base_plant.leaf.k_leaf_max
    for i in tqdm(list(range(len(f_range))), desc="param_it"):
        plant = base_plant
        k_actual_trial = np.zeros(num_s)
        k_actual_frac_trial = np.zeros(num_s)
        k_stem_frac_trial = np.zeros(num_s)
        k_stem_trial = np.zeros(num_s)
        counter = 0
        plant.stem.k_stem_max = k_stem_max * plant.stem.f_stem / f_range[i]
        plant.leaf.k_leaf_max = k_leaf_max * plant.leaf.f_leaf / (1 - f_range[i])
        for j in range(len(s_range)):
            s0 = s_range[j]
            k_stem, k_actual, k_stem_frac, k_actual_frac, _, _ = ss_simulate_psi_k_norefilling(s0, plant, VPD)
            k_actual_trial[counter] = k_actual
            k_stem_trial[counter] = k_stem
            k_stem_frac_trial[counter] = k_stem_frac
            k_actual_frac_trial[counter] = k_actual_frac
            counter = counter + 1

        np.savetxt(os.path.join(path, "k_plant_mean.txt"), k_actual_trial, delimiter=',')
        np.savetxt(os.path.join(path, "k_stem_mean.txt"), k_stem_trial, delimiter=',')
        np.savetxt(os.path.join(path, "k_plant_frac_mean.txt"), k_actual_frac_trial, delimiter=',')
        np.savetxt(os.path.join(path, "k_stem_frac_mean.txt"), k_stem_frac_trial, delimiter=',')
        np.savetxt(os.path.join(path, "psi_s.txt"), s_range, delimiter=',')


if __name__ == '__main__':
    plot_ss_soil_vs_k(pop, VPD, s_range)


