import numpy as np
import scipy.optimize as opt
import params_constants


# function for clipping an input value between designated minimum and maximum
clamp = lambda minimum, x, maximum: max(minimum, min(x, maximum))

# Threshold parameters

pxf = 5  # a multiplier to P50 beyond which xylem water potential is not allowed
plf = 5  # a multiplier to P50 beyond which leaf water potential is not allowed


class Canopy:
    def __init__(self, Gs_leaf, A_canopy, Gs_min, H_canopy):
        self.Gs_leaf = Gs_leaf  # max leaf conductance for H2O, mol/m2/s
        self.A_canopy = A_canopy  # canopy area, m2
        self.H_canopy = H_canopy  # Mean canopy height, m
        self.Gs_min = Gs_min  # min leaf conductance for H2O, mol/m2/s

    def g_canopy(self, gs):
        # Maximum canopy conductance per day, m3/d
        Sd, Vw = params_constants.Sd, params_constants.Vw
        return gs * self.A_canopy * Vw * Sd


class Stem:
    def __init__(self, P50_stem, a_stem, k_plant_max, f_stem):
        self.P50_stem = P50_stem  # water potential at 50% loss in maximum xylem conductivity, MPa
        self.a_stem = a_stem  # parameter of xylem vulnerability in 'plc' form
        self.f_stem = f_stem
        self.k_stem_max = k_plant_max / f_stem  # fraction of plant conductance in stem

    def k_stem(self, P_stem):
        # Stem conductance, m3/d-MPa
        kmax, a, P50 = self.k_stem_max, self.a_stem, self.P50_stem
        # based on exponential-sigmoidal function
        # constrain the domain of P to be negative
        return kmax * (1.0 - 1.0 / (1.0 + np.exp(a * (clamp(pxf * P50, P_stem, 0) - P50))))

    def k_stem_Kirchhoff(self, P):
        # Kirchoff transform - integrate vulnerability curve from -infinity to P_stem (P)
        kmax, a, P50 = self.k_stem_max, self.a_stem, self.P50_stem
        # based on exponential-sigmoidal function
        # constrain the domain of P to be negative
        k_stem = kmax * (np.log(np.exp(-a * P50) + np.exp(- a * clamp(pxf * P50, P, 0))) / a + P)
        return k_stem

    def Px_solvent(self, P_stem, percent):
        # helper function for solving for P_stem corresponding to designated percent loss of conductivity
        return self.k_stem(P_stem) - (1.0 - percent) * self.k_stem_max


class Leaf:
    def __init__(self, P50_leaf, a_leaf, k_max_total, f_leaf):
        self.P50_leaf = P50_leaf  # water potential at 50% loss in maximum leaf conductivity, MPa
        self.a_leaf = a_leaf  # parameter of leaf vulnerability in 'plc' form
        self.f_leaf = f_leaf
        self.k_leaf_max = k_max_total / f_leaf  # fraction of plant conductance in stem

    def k_leaf(self, P_leaf1):
        # Leaf conductance, m3/d-MPa
        kmax, a, P50 = self.k_leaf_max, self.a_leaf, self.P50_leaf
        # based on exponential-sigmoidal function
        # constrain the domain of P to be negative
        return kmax * (1.0 - 1.0 / (1.0 + np.exp(a * (clamp(plf * P50, P_leaf1, 0) - P50))))

    def k_leaf_Kirchhoff(self, P):
        # Kirchoff transform - integrate vulnerability curve from -infinity to P_leaf (P)
        kmax, a, P50 = self.k_leaf_max, self.a_leaf, self.P50_leaf
        # based on exponential-sigmoidal function
        # constrain the domain of P to be negative
        k_leaf = kmax * (np.log(np.exp(-a * P50) + np.exp(- a * clamp(plf * P50, P, 0))) / a + P)
        return k_leaf

    def Pl_solvent(self, P_leaf1, percent):
        # helper function for solving for P_stem corresponding to designated percent loss of conductivity
        return self.k_leaf(P_leaf1) - (1.0 - percent) * self.k_leaf_max


class Whole_plant:
    def __init__(self, species):
        self.species = species  # Species label
        self.Vw = 18 * 10 ** (-6)  # Partial molal volume of water m3/mol

    @property
    def canopy(self):
        return self._canopy

    @canopy.setter
    def canopy(self, canopy_object):
        self._canopy = canopy_object

    @property
    def stem(self):
        return self._stem

    @stem.setter
    def stem(self, stem_object):
        self._stem = stem_object

    @property
    def roots(self):
        return self._roots

    @roots.setter
    def roots(self, roots_object):
        self._roots = roots_object

    @property
    def leaf(self):
        return self._leaf

    @leaf.setter
    def leaf(self, leaf_object):
        self._leaf = leaf_object

    def stomata_func(self, P_leaf2):
        # Stomatal conductance as a funtion of leaf water potential
        gs_max, gs_min, P50_stem = self.canopy.Gs_leaf, self.canopy.Gs_min, self.stem.P50_stem
        # find Pl_90 based on relationship between P50_stem and Pclose from Martin St Paul et al (2017)
        if P50_stem >= -6.69:
            pl_90 = 0.3231 * P50_stem - 0.9264
        else:
            pl_90 = 0.0052 * P50_stem - 3.0516
        gs1 = gs_max * (1 - .9 * (P_leaf2 / pl_90))
        gs = max(gs_min, gs1)
        return gs

    def ET_soil_func_fractional(self, P_soil, P_stem, P_leaf1):
        # calcualtes soil-root conductance as constant fraction of resistance from soil to leaf (see Tyree and Sperry (1988))
        k_stem, k_leaf = self.stem.k_stem, self.leaf.k_leaf
        frac_resistance = 0.18
        if k_stem(P_stem) >= 0.00001 and k_leaf(P_leaf1) >= 0.00001:
            k_sr = ((1.0 - frac_resistance) / frac_resistance) / ((1.0 / k_stem(P_stem)) + (1.0 / k_leaf(P_leaf1)))
        elif k_stem(P_stem) < 0.00001 and k_leaf(P_leaf1) >= 0.00001:
            k_sr = ((1.0 - frac_resistance) / frac_resistance) / ((1.0 / 0.00001) + (1.0 / k_leaf(P_leaf1)))
        elif k_stem(P_stem) >= 0.00001 and k_leaf(P_leaf1) < 0.00001:
            k_sr = ((1.0 - frac_resistance) / frac_resistance) / ((1.0 / k_stem(P_stem)) + (1.0 / 0.00001))
        else:
            k_sr = ((1.0 - frac_resistance) / frac_resistance) / ((1.0 / 0.00001) + (1.0 / 0.00001))
        return k_sr * (P_soil - P_stem)

    def estimate_fromPx_ET_soil_func_fractional(self, P_soil, P_stem):
        # calcualtes soil-root conductance as constant fraction of resistance from soil to leaf (see Tyree and Sperry (1988))
        # to estimate for initial guess, assume kstem=kleaf
        k_stem = self.stem.k_stem
        frac_resistance = 0.18
        if k_stem(P_stem) > 0.00001:
            k_sr = ((1.0 - frac_resistance) / frac_resistance) / (
                        (1.0 / k_stem(P_stem)) + (1.0 / k_stem(k_stem(P_stem))))
        else:
            k_sr = 0.0
        return k_sr * (P_soil - P_stem)

    def ET_stem_func(self, P_stem, P_leaf, Px_min):
        # Transpiration from stem to canopy, m3/d; Equation (1) of Feng (2018)
        k_stemKirch = self.stem.k_stem_Kirchhoff
        k_stem, P50 = self.stem.k_stem, self.stem.P50_stem
        Px_clamped = clamp(pxf * P50, P_stem, 0)
        # To account for refilling
        return k_stemKirch(Px_clamped) - k_stemKirch(min(P_leaf, P_stem)) + max(0, k_stem(Px_min) * (P_stem - Px_min))

    def ET_leaf_func(self, P_leaf1, P_leaf2, Pl_min):
        # Transpiration through leaf tissues, m3/d; Equation (1) of Feng (2018)
        k_leafKirch = self.leaf.k_leaf_Kirchhoff
        k_leaf, P50 = self.leaf.k_leaf, self.leaf.P50_leaf
        Pl_clamped = clamp(plf * P50, P_leaf1, 0)
        # To account for refilling
        return k_leafKirch(Pl_clamped) - k_leafKirch(min(P_leaf2, P_leaf1)) + max(0,
                                                                                  k_leaf(Pl_min) * (P_leaf1 - Pl_min))

    def ET_canopy_func(self, gs, VPD):
        # Transpiration from the canopy; m3/d
        g_canopy = self.canopy.g_canopy
        P_baro = params_constants.P_baro
        return g_canopy(gs) * (VPD / P_baro)

    def Pl1_ET_func(self, ET, P_stem, Px_min):
        # Solving for leaf water potential Pl given (ET, P_stem, and Px_min)
        # Px_min is the minimum stem water potential experienced by the plant
        a, P50, kmax, k_stem = self.stem.a_stem, self.stem.P50_stem, self.stem.k_stem_max, self.stem.k_stem

        # change the effective ET and P_stem processed over the untraversed area of the vulnerability curve (v.c.) to account for no refilling
        ET_eff = ET - max(0, k_stem(Px_min) * (
                    P_stem - Px_min))  # effective ET; if P_stem > Px_min, then account for area over kmax
        # if ET_eff< 0:
        # ET_eff=ET
        P_stem = min(P_stem, Px_min)  # effective P_stem; if P_stem > Px_min, then start traversing v.c. from Px_min

        # vulnerability curve in "plc" form
        h = 1.0 - np.exp(a * ET_eff / kmax) + np.exp(a * (P_stem - P50))
        if h <= 0:
            Pl1 = 1.1 * P_stem
        else:
            Pl1 = -ET_eff / kmax + P50 + 1.0 / a * np.log(h)
        return Pl1

    def Pl2_ET_func(self, ET, P_leaf1, Pl_min):
        # Solving for leaf water potential Pl2 given (ET, P_leaf1, and Pl_min)
        # Px_min is the minimum stem water potential experienced by the plant
        a, P50, kmax, k_leaf = self.leaf.a_leaf, self.leaf.P50_leaf, self.leaf.k_leaf_max, self.leaf.k_leaf

        # change the effective ET and P_leaf1 processed over the untraversed area of the vulnerability curve (v.c.) to account for no refilling
        ET_eff = ET - max(0, k_leaf(Pl_min) * (
                    P_leaf1 - Pl_min))  # effective ET; if P_leaf1 > Pl_min, then account for area over kmax
        # if ET_eff< 0:
        # ET_eff=ET
        P_leaf1 = min(P_leaf1, Pl_min)  # effective P_stem; if P_leaf1 > Pl_min, then start traversing v.c. from Pl_min

        # vulnerability curve in "plc" form
        h = 1.0 - np.exp(a * ET_eff / kmax) + np.exp(a * (P_leaf1 - P50))
        if h <= 0:
            Pl2 = 1.1 * P_leaf1
        else:
            Pl2 = -ET_eff / kmax + P50 + 1.0 / a * np.log(h)
        return Pl2

    def flux_solvent(self, xxx_todo_changeme, P_soil, VPD, Px_min, Pl_min):
        # Helper function for solving coupled soil-plant hydraulic system
        # A can be solved later as function of gs - Equation (9) of Feng (2018)
        # A - gs*kcarbx *CO2conc/(1.6*kcarbx + gs)     # mol/m2-s
        # Constraints: gs need to be greatern than 0, minimum P_stem
        (ET, P_stem, P_leaf1, gs, P_leaf2) = xxx_todo_changeme
        gs_min = self.canopy.Gs_min

        return (max(gs, gs_min) - self.stomata_func(P_leaf2),
                ET - self.ET_soil_func_fractional(P_soil, P_stem, P_leaf1),
                ET - abs(self.ET_stem_func(P_stem, P_leaf1, Px_min)),
                ET - abs(self.ET_leaf_func(P_leaf1, P_leaf2, Pl_min)),
                ET - abs(self.ET_canopy_func(gs, VPD)))

    def flux_solver(self, P_soil, VPD, Px_min, Pl_min, init):
        # Solves coupled soil-plant hydraulic system - system of equations describing internal water states (ET, P_stem, P_leaf1, P_leaf2)
        # based on an input of soil water potential and VPD
        sol = opt.root(self.flux_solvent, init, (P_soil, VPD, Px_min, Pl_min), method='hybr')
        return sol

    def calculate_from_Px(self, Px, VPD, Px_min, Ps, Pl_min):
        # Calculates system variables based on designated stem water potential Px.
        J = self.estimate_fromPx_ET_soil_func_fractional(Ps, Px)  # Transpiration from soil to stem
        Pl1 = self.Pl1_ET_func(J, Px, Px_min)  # Leaf water potential 1
        Pl2 = self.Pl2_ET_func(J, Pl1, Pl_min)  # Leaf water potential 2
        gs = self.stomata_func(Pl2)  # Stomatal conductance
        ET = self.ET_canopy_func(gs, VPD)  # Transpiration from canopy
        return J, Pl1, Pl2, gs, ET

    def get_init_conditions_psi(self, VPD, Ps, Px_min, Pl_min, Px_arr):
        # Outputs guesses for initial conditions used to solve the system, based on inputs of soil water potential, VPD, and minimum stem water potential Px_min
        # Px_arr is the range of stem water potentials over which potential initial conditions are tested.
        Inits_arr = np.zeros((len(Px_arr), 5));
        Inits_arr2 = np.zeros_like(Inits_arr)
        diffET_arr = np.zeros_like(Px_arr);
        J_arr = np.zeros_like(Px_arr);
        ET_arr = np.zeros_like(Px_arr)
        for i, Px in enumerate(Px_arr):
            J, Pl1, Pl2, gs, ET = self.calculate_from_Px(Px, VPD, Px_min, Ps, Pl_min)
            Inits_arr[i] = [ET, Px, Pl1, gs, Pl2]  # were these variables out of order?
            Inits_arr2[i] = [J, Px, Pl1, gs, Pl2]
            J_arr[i] = J
            ET_arr[i] = ET
            diffET_arr[i] = J - ET
        # find when difference between water supply and demand is minimal
        j = np.argmin(np.abs(diffET_arr))
        init_list = []
        try:
            for k in [0, 1, -1]:
                init_list.extend([Inits_arr[j + k]])
                init_list.extend([Inits_arr2[j + k]])
        except IndexError:  # if minimum difference for water supply and demand is found close to the edge of Px range
            init_list.extend([Inits_arr[-1], Inits_arr2[-1]])
            init_list.extend([Inits_arr[0], Inits_arr2[0]])
        return init_list


    def get_fluxes_scalar_psi(self, VPD, P_soil, Px_min, Pl_min, Flux_inits=None):
        # returns solutions to the soil-plant system (E, P_stem, P_leaf1, P_leaf2)

        # Set a range of stem water potentials Px over which initial conditions can be guessed
        Px_arr = np.linspace(P_soil * 1.2, P_soil, 50)

        sol_flag = False  # flag will be turned to True only when solution is found with inputted guesses
        if Flux_inits is not None:  # if initial conditions have already been included in the arguments
            sol0 = self.flux_solver(P_soil, VPD, Px_min, Pl_min, Flux_inits)
            sol_flag = sol0.success
            if sol_flag:
                ET, P_stem, P_leaf1, gs, P_leaf2 = sol0.x
        # if initial conditions have not been designated, or if previous attempt at solving was not successful,
        # then we used guessed initial conditions to solve the system
        if (Flux_inits is None) or (sol_flag is False):
            Flux_inits_generated = self.get_init_conditions_psi(VPD, P_soil, Px_min, Pl_min, Px_arr)
            sol1 = self.flux_solver(P_soil, VPD, Px_min, Pl_min, Flux_inits_generated[0])

            if sol1.success:
                ET, P_stem, P_leaf1, gs, P_leaf2 = sol1.x
            else:
                for _, init in enumerate(Flux_inits_generated[1:]):
                    sol = self.flux_solver(P_soil, VPD, Px_min, Pl_min, init)
                    if sol.success:
                        ET, P_stem, P_leaf1, gs, P_leaf2 = sol.x
                        break
                    # the default is the best initial guess
                    ET, P_stem, P_leaf1, gs, P_leaf2 = Flux_inits_generated[0]
        # calculate assimilation independently
        # return variables
        return [ET, P_stem, P_leaf1, P_leaf2, gs]
