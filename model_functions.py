import xlrd
import numpy as np
import scipy.stats as st
import pandas as pd
import params_constants
from SPAC_model import Canopy, Stem, Whole_plant, Leaf


# adjust this one after making traits file
def import_traits_data(filepath=None, species=None, sp_coln=None):
    if filepath == None: filepath = './hydraulic_traits.xls'
    if species == None: species = ['POP', 'OAK', 'TEST', 'GENERIC', 'POP_REV', 'OAK_REV']
    if sp_coln == None: sp_coln = [2, 4, 6, 8, 10, 12]
    book = xlrd.open_workbook(filepath)
    sheet = book.sheet_by_name('parameters')
    keys = np.asarray([_f for _f in sheet.col_values(0) if _f], dtype='str')
    canopy_keys = ['A_canopy', 'Gs_leaf', 'Gs_min', 'H_canopy']
    stem_keys = ['f_stem', 'k_plant_max', 'a_stem', 'P50_stem']
    root_keys = ['L_root', 'A_root', 'd_root']
    leaf_keys = ['a_leaf', 'P50_leaf', 'f_leaf', 'k_max_total']
    chap_dict = {}
    for sp, spc in zip(species, sp_coln):
        chap_dict[sp] = {}
        for part_keys, part_dict in zip([canopy_keys, stem_keys, root_keys, leaf_keys],
                                        ['canopy_dict', 'stem_dict', 'root_dict', 'leaf_dict']):
            chap_dict[sp][part_dict] = {}
            for key in part_keys:
                j = np.where(keys == key)[0][0] + 1  # to account for the column of the species
                chap_dict[sp][part_dict][key] = sheet.col_values(spc)[j]
    return chap_dict


def initialize_plant(sp, params, soil_type=None):
    # initializing each species with its name (in strings)
    canopy_dict = params[sp]['canopy_dict']
    stem_dict = params[sp]['stem_dict']
    root_dict = params[sp]['root_dict']
    leaf_dict = params[sp]['leaf_dict']
    plant = Whole_plant(species=sp)
    plant.canopy = Canopy(**canopy_dict)
    plant.stem = Stem(**stem_dict)
    plant.leaf = Leaf(**leaf_dict)
    # plant.roots = Roots(soil_type=soil_type, **root_dict)
    return plant


def get_part(var):
    if var in ['A_canopy', 'Gs_leaf', 'Amax', 'rho', 'lamABA', 'm']:
        return 'canopy_dict'
    elif var in ['L_stem', 'A_stem', 'Ksat_stem', 'a_stem', 'c_stem', 'P50_stem']:
        return 'stem_dict'
    elif var in ['L_root', 'A_root']:
        return 'root_dict'


def simulate_s_t_norefilling(tRun, dt, sInit, plant, VPD):
    s_t = np.zeros(len(tRun));
    s0 = sInit
    px_t = np.zeros_like(s_t)
    pxmin_t = np.zeros_like(s_t)
    E_t = np.zeros_like(s_t)
    gs_t = np.zeros_like(s_t)
    pl1_t = np.zeros_like(s_t)
    pl2_t = np.zeros_like(s_t)
    plmin_t = np.zeros_like(s_t)

    Px_min = plant.get_fluxes_scalar(VPD, s0, 0, 0)[1]
    Pl_min = Px_min
    px_t[0] = Px_min

    for i in range(len(tRun)):

        # Hold soil moisture constant
        s1 = s0
        P_soil = plant.roots.P_soil_solver(s1)

        # Calculate fluxes and plant water potentials
        if i > 0:
            Px = plant.get_fluxes_scalar(VPD, s0, Px_min, Pl_min, prev_sol)[1]
        else:
            Px = plant.get_fluxes_scalar(VPD, s0, Px_min, Pl_min)[1]
        J = plant.ET_soil_func(P_soil, Px)
        Pl1 = plant.Pl1_ET_func(J, Px, Px_min)
        Pl2 = plant.Pl2_ET_func(J, Pl1, Pl_min)

        # calculate the rest of the variables
        gs = plant.stomata_func(Pl2)
        E = plant.ET_canopy_func(gs, VPD)

        # check Px and Pl to see if it has exceeded Px_min or Pl_min
        # Px_min = min(Px, Px_min)
        # Pl_min = min(Pl1, Pl_min)

        Px_min = Px
        Pl_min = Pl1

        # use previous step as initial conditions for next step
        prev_sol = [E, Px, Pl1, gs, Pl2]

        # update to next step
        s_t[i] = s0;
        px_t[i] = Px
        pxmin_t[i] = Px_min
        E_t[i] = E
        gs_t[i] = gs
        pl1_t[i] = Pl1
        pl2_t[i] = Pl2
        plmin_t[i] = Pl_min
        # return variables
    return s_t, px_t, pxmin_t, E_t, gs_t, pl1_t, pl2_t, plmin_t


def ss_simulate_s_t_norefilling(s0, plant, VPD):
    max_iterations = 200
    tolerance = .001
    converged = 0
    min_iterations = 5

    s_stepup = s0 * 1.05
    Px_min = plant.get_fluxes_scalar(VPD, s_stepup, 0, 0)[1]
    Pl_min = Px_min
    # Hold soil moisture constant
    P_soil = plant.roots.P_soil_solver(s0)

    for i in range(max_iterations):

        # Calculate fluxes and plant water potentials
        if i > 0:
            Px = plant.get_fluxes_scalar(VPD, s0, Px_min, Pl_min, prev_sol)[1]
        else:
            Px = plant.get_fluxes_scalar(VPD, s0, Px_min, Pl_min)[1]

        # checks to see if Px value has converged to within tolerance
        if i > 0:
            Px_change = abs(Px - Px_prev) / abs(Px_prev)
            if Px_change <= tolerance:
                converged = 1
            else:
                converged = 0

        J = plant.ET_soil_func(P_soil, Px)
        Pl1 = plant.Pl1_ET_func(J, Px, Px_min)
        Pl2 = plant.Pl2_ET_func(J, Pl1, Pl_min)

        # calculate the rest of the variables
        gs = plant.stomata_func(Pl2)
        E = plant.ET_canopy_func(gs, VPD)

        # check Px and Pl to see if it has exceeded Px_min or Pl_min
        # Px_min = min(Px, Px_min)
        # Pl_min = min(Pl1, Pl_min)

        Px_min = (Px + 2.0 * Px_min) / 3.0
        Pl_min = (Pl1 + 2.0 * Pl_min) / 3.0

        # use previous step as initial conditions for next step
        prev_sol = [E, Px, Pl1, gs, Pl2]

        Px_prev = Px

        if converged == 1 & i >= min_iterations:  # breaks loop early if it has reached convergence, otherwise will loop to max iterations
            break
        if i == max_iterations:
            print(("VPD =", VPD, "s=", s0, "Px_change=", Px_change))

    k_stem = plant.stem.k_stem(Px)
    k_leaf = plant.leaf.k_leaf(Pl1)
    k_actual = (k_stem * k_leaf) / (k_stem + k_leaf)
    original_k_max = (plant.stem.k_stem_max * plant.leaf.k_leaf_max) / (plant.stem.k_stem_max + plant.leaf.k_leaf_max)
    k_actual_frac = k_actual / original_k_max
    k_tissue_change = (k_leaf / plant.leaf.k_leaf_max) / (k_stem / plant.stem.k_stem_max)

    # return variables
    return k_stem, k_actual, k_actual_frac, k_tissue_change, Px


def ss_simulate_psi_k_norefilling(psi_s0, plant, VPD):
    max_iterations = 200
    tolerance = .001
    converged = 0
    min_iterations = 5

    psi_s_stepup = psi_s0 * 0.9
    Px_min = plant.get_fluxes_scalar_psi(VPD, psi_s0, psi_s0, psi_s0)[1]
    Pl_min = Px_min

    for i in range(max_iterations):

        # Calculate fluxes and plant water potentials
        if i > 0:
            _, Px, Pl1, _, _ = plant.get_fluxes_scalar_psi(VPD, psi_s0, Px_min, Pl_min, prev_sol)
        else:
            _, Px, Pl1, _, _ = plant.get_fluxes_scalar_psi(VPD, psi_s0, Px_min, Pl_min)

        # checks to see if Px value has converged to within tolerance
        if i > 0:
            Px_change = abs(Px - Px_prev) / abs(Px_prev)
            if Px_change <= tolerance:
                converged = 1
            else:
                converged = 0

        J = plant.ET_soil_func_fractional(psi_s0, Px, Pl1)
        Pl1 = plant.Pl1_ET_func(J, Px, Px_min)
        Pl2 = plant.Pl2_ET_func(J, Pl1, Pl_min)

        # calculate the rest of the variables
        gs = plant.stomata_func(Pl2)
        E = plant.ET_canopy_func(gs, VPD)

        # check Px and Pl to see if it has exceeded Px_min or Pl_min
        # Px_min = min(Px, Px_min)
        # Pl_min = min(Pl1, Pl_min)

        Px_min = (Px + 2.0 * Px_min) / 3.0
        Pl_min = (Pl1 + 2.0 * Pl_min) / 3.0

        # use previous step as initial conditions for next step
        prev_sol = [E, Px, Pl1, gs, Pl2]

        Px_prev = Px

        if converged == 1 & i >= min_iterations:  # breaks loop early if it has reached convergence, otherwise will loop to max iterations
            break
        if i == max_iterations:
            print(("VPD =", VPD, "psi_s=", psi_s0, "Px_change=", Px_change))

    k_stem = plant.stem.k_stem(Px)
    k_leaf = plant.leaf.k_leaf(Pl1)
    k_actual = (k_stem * k_leaf) / (k_stem + k_leaf)
    original_k_max = (plant.stem.k_stem_max * plant.leaf.k_leaf_max) / (plant.stem.k_stem_max + plant.leaf.k_leaf_max)
    k_stem_frac = k_stem / plant.stem.k_stem_max
    k_actual_frac = k_actual / original_k_max
    k_tissue_change = (k_leaf / plant.leaf.k_leaf_max) / (k_stem / plant.stem.k_stem_max)

    # return variables
    return k_stem, k_actual, k_stem_frac, k_actual_frac, k_tissue_change, Px


def ss_simulate_psi_k_norefilling_test(psi_s0, plant, VPD):
    max_iterations = 200
    tolerance = .001
    converged = 0
    min_iterations = 5

    psi_s_stepup = psi_s0 * 0.9
    Px_min = plant.get_fluxes_scalar_psi(VPD, psi_s0, psi_s0, psi_s0)[1]
    Pl_min = Px_min

    for i in range(max_iterations):

        # Calculate fluxes and plant water potentials
        if i > 0:
            _, Px, Pl1, _, _ = plant.get_fluxes_scalar_psi(VPD, psi_s0, Px_min, Pl_min, prev_sol)
        else:
            _, Px, Pl1, _, _ = plant.get_fluxes_scalar_psi(VPD, psi_s0, Px_min, Pl_min)

        # checks to see if Px value has converged to within tolerance
        if i > 0:
            Px_change = abs(Px - Px_prev) / abs(Px_prev)
            if Px_change <= tolerance:
                converged = 1
            else:
                converged = 0

        J = plant.ET_soil_func_fractional(psi_s0, Px, Pl1)
        Pl1 = plant.Pl1_ET_func(J, Px, Px_min)
        Pl2 = plant.Pl2_ET_func(J, Pl1, Pl_min)

        # calculate the rest of the variables
        gs = plant.stomata_func(Pl2)
        E = plant.ET_canopy_func(gs, VPD)

        # check Px and Pl to see if it has exceeded Px_min or Pl_min
        # Px_min = min(Px, Px_min)
        # Pl_min = min(Pl1, Pl_min)

        Px_min = (Px + 2.0 * Px_min) / 3.0
        Pl_min = (Pl1 + 2.0 * Pl_min) / 3.0

        # use previous step as initial conditions for next step
        prev_sol = [E, Px, Pl1, gs, Pl2]

        Px_prev = Px

        if converged == 1 & i >= min_iterations:  # breaks loop early if it has reached convergence, otherwise will loop to max iterations
            break
        if i == max_iterations:
            print(("VPD =", VPD, "psi_s=", psi_s0, "Px_change=", Px_change))

    k_stem_1 = plant.stem.k_stem(Px)
    k_stem_2 = plant.stem.k_stem(Pl1)
    k_stem_3 = plant.stem.k_stem(((Pl1 + Px) / 2.0))
    k_stem_4 = E / (Px - Pl1)
    k_leaf = plant.leaf.k_leaf(Pl1)
    k_stem_new = E / (Px - Pl1)
    k_leaf_new = E / (Pl1 - Pl2)
    k_actual = (k_stem_1 * k_leaf) / (k_stem_1 + k_leaf)
    k_actual_new = (k_stem_new * k_leaf_new) / (k_stem_new + k_leaf_new)
    original_k_max = (plant.stem.k_stem_max * plant.leaf.k_leaf_max) / (plant.stem.k_stem_max + plant.leaf.k_leaf_max)
    k_stem_frac = k_stem_1 / plant.stem.k_stem_max
    k_actual_frac = k_actual / original_k_max
    k_tissue_change = (k_leaf / plant.leaf.k_leaf_max) / (k_stem_1 / plant.stem.k_stem_max)

    # return variables
    return k_stem_1, k_stem_2, k_stem_3, k_stem_4, k_actual, k_stem_new, k_actual_new, Px, Pl1, Pl2


def gen_k_from_P50(P50_range, Leaf_Area, Lx, num_sims):
    XFT = pd.read_csv('cleaned_XFT.txt', sep='\t')
    XFT_Huber = pd.read_csv('cleaned_XFT_Huber.txt', sep='\t')
    XFT_slope = pd.read_csv('cleaned_XFT_slope.txt', sep='\t')
    min_P50 = np.min(XFT.P50)
    max_P50 = np.max(XFT.P50)
    Sd, rhoH2O = params_constants.Sd, params_constants.rhoH2O
    k_min = np.min(XFT.k_sat)
    k_max = np.max(XFT.k_sat)
    k_vals_range = np.linspace(k_min, k_max, 500)
    data = np.vstack([XFT.P50, XFT.k_sat])
    kde = st.gaussian_kde(data)
    slope_min = np.min(XFT_slope.Slope)
    slope_max = np.max(XFT_slope.Slope)
    slope_vals_range = np.linspace(slope_min, slope_max, 500)
    data_slope = np.vstack([XFT_slope.P50, XFT_slope.Slope])
    kde_slope = st.gaussian_kde(data_slope)
    for i in range(len(P50_range)):
        P50 = -P50_range[i]  # need to convert to positive value
        P50_vals = np.linspace(P50, P50, 500)
        Z = kde.evaluate(np.vstack([P50_vals, k_vals_range]))
        Z_norm = Z / np.sum(Z)
        Z_counts = 100000.0 * Z_norm
        Z_counts = Z_counts.astype(int)
        k_vals = np.repeat(k_vals_range, Z_counts)
        params_dist = st.lognorm.fit(k_vals, floc=0)
        P50_simk = st.lognorm.rvs(params_dist[0], loc=params_dist[1], scale=params_dist[2], size=num_sims)
        Z_slope = kde_slope.evaluate(np.vstack([P50_vals, slope_vals_range]))
        Z_slope_norm = Z_slope / np.sum(Z_slope)
        Z_slope_counts = 100000.0 * Z_slope_norm
        Z_slope_counts = Z_slope_counts.astype(int)
        slope_vals = np.repeat(slope_vals_range, Z_slope_counts)
        # params_dist = st.lognorm.fit(k_vals, floc=0)
        # P50_simk = st.lognorm.rvs(params_dist[0], loc=params_dist[1], scale=params_dist[2], size=num_sims)
        params_slope_dist = st.lognorm.fit(slope_vals, floc=0)
        P50_simslope = st.lognorm.rvs(params_slope_dist[0], loc=params_slope_dist[1], scale=params_slope_dist[2],
                                      size=num_sims)
        P50_sim_a = -np.log((1.0 / .88) - 1) / P50_simslope * (88.0 - 50.0)  # convert slope to a parameter for VC
        if i == 0:
            all_simk = P50_simk
            all_sima = P50_sim_a
            all_P50 = np.linspace(P50, P50, num_sims)
        else:
            all_simk = np.append(all_simk, P50_simk)
            all_sima = np.append(all_sima, P50_sim_a)
            new_P50 = np.linspace(P50, P50, num_sims)
            all_P50 = np.append(all_P50, new_P50)
    H_data = np.vstack([XFT_Huber.k_sat, XFT_Huber.Huber])
    H_kde = st.gaussian_kde(H_data)
    H_min = np.min(XFT_Huber.Huber)
    H_max = np.max(XFT_Huber.Huber)
    H_vals_range = np.linspace(H_min, H_max, 500)
    for i in range(len(all_simk)):
        k = all_simk[i]
        k_orig = k
        if k > 12.0:
            k = 12.0
        k_vals = np.linspace(k, k, 500)
        Z = H_kde.evaluate(np.vstack([k_vals, H_vals_range]))
        Z_norm = Z / np.sum(Z)
        Z_counts = 100000.0 * Z_norm
        Z_counts = Z_counts.astype(int)
        H_vals = np.repeat(H_vals_range, Z_counts)
        params_dist = st.gamma.fit(H_vals, floc=0)
        k_simH = st.gamma.rvs(params_dist[0], loc=params_dist[1], scale=params_dist[2], size=1)
        sim_cond = k_orig * k_simH * Leaf_Area * Sd / (Lx * rhoH2O)
        if i == 0:
            all_simH = k_simH
            all_simcond = sim_cond
        else:
            all_simH = np.append(all_simH, k_simH)
            all_simcond = np.append(all_simcond, sim_cond)
    return all_simcond, all_sima