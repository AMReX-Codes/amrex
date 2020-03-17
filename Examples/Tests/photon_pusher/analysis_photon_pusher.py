#! /usr/bin/env python

# Copyright 2019 Luca Fedeli, Maxence Thevenet, Weiqun Zhang
#
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

import yt
import numpy as np
import sys

#This script checks if photons initialized with different momenta and
#different initial directions propagate along straight lines at the speed of
#light. The plotfile to be analyzed is passed as a command line argument.

#If the script is run without a command line argument, it regenerates a new
#inputfile according to the initial conditions listed below.


#Physical constants
c = 299792458.
m_e = 9.1093837015e-31
#________________________________________

#Test cases
spec_names = ["p_xp_1", "p_xn_1", "p_yp_1", "p_yn_1",
    "p_zp_1", "p_zn_1","p_dp_1", "p_dn_1",
    "p_xp_10", "p_xn_10", "p_yp_10", "p_yn_10",
    "p_zp_10", "p_zn_10", "p_dp_10", "p_dn_10"]
#photon momenta are in units of m_e c
mxp1 = np.array([1, 0.0, 0.0])
mxn1 = np.array([-1, 0.0, 0.0])
myp1 = np.array([0.0, 1, 0.0])
myn1 = np.array([0.0, -1, 0.0])
mzp1 = np.array([0.0, 0.0, 1])
mzn1 = np.array([0.0, 0.0, -1])
mdp1 = np.array([1, 1, 1])
mdn1 = np.array([-1, -1, -1])
mxp10 = np.array([10, 0.0, 0.0])
mxn10 = np.array([-10, 0.0, 0.0])
myp10 = np.array([0.0, 10, 0.0])
myn10 = np.array([0.0, -10, 0.0])
mzp10 = np.array([0.0, 0.0, 10])
mzn10 = np.array([0.0, 0.0, -10])
mdp10 = np.array([10, 10, 10])
mdn10 = np.array([-10,-10, -10])
gamma_beta_list = np.array([mxp1, mxn1, myp1, myn1, mzp1, mzn1, mdp1, mdn1,
    mxp10, mxn10, myp10, myn10, mzp10, mzn10, mdp10, mdn10])
init_pos =  np.array([0.0, 0.0, 0.0])
#________________________________________

#Tolerance
tol_pos = 1.0e-14;
tol_mom = 0.0; #momentum should be conserved exactly
#________________________________________

#Input filename
inputname = "inputs"
#________________________________________

# This function reads the WarpX plotfile given as the first command-line
# argument, and check if the position of each photon agrees with theory.
def check():
    filename = sys.argv[1]
    data_set_end = yt.load(filename)

    sim_time = data_set_end.current_time.to_value()

    #expected positions list
    ll = sim_time*c
    answ_pos = init_pos + \
    ll*gamma_beta_list/np.linalg.norm(gamma_beta_list,axis=1, keepdims=True)

    #expected momenta list
    answ_mom =  m_e * c *gamma_beta_list #momenta don't change

    #simulation results
    all_data =  data_set_end.all_data()
    res_pos = [np.array([
        all_data[sp, 'particle_position_x'].v[0],
        all_data[sp, 'particle_position_y'].v[0],
        all_data[sp, 'particle_position_z'].v[0]])
        for sp in spec_names]
    res_mom = [np.array([
        all_data[sp, 'particle_momentum_x'].v[0],
        all_data[sp, 'particle_momentum_y'].v[0],
        all_data[sp, 'particle_momentum_z'].v[0]])
        for sp in spec_names]

    #check discrepancies
    disc_pos = [np.linalg.norm(a-b)/np.linalg.norm(b)
        for a,b in zip(res_pos, answ_pos)]
    disc_mom = [np.linalg.norm(a-b)/np.linalg.norm(b)
        for a,b in zip(res_mom, answ_mom)]

    print("max(disc_pos) = %s" %max(disc_pos))
    print("tol_pos = %s" %tol_pos)
    print("max(disc_mom) = %s" %max(disc_mom))
    print("tol_mom = %s" %tol_mom)

    assert ((max(disc_pos) <= tol_pos) and (max(disc_mom) <= tol_mom))

# This function generates the input file to test the photon pusher.
def generate():
    with open(inputname,'w') as f:
        f.write("#Automatically generated inputfile\n")
        f.write("#Run check.py without arguments to regenerate\n")
        f.write("#\n\n")
        f.write("max_step = 50\n")
        f.write("amr.n_cell = 64 64 64\n")
        f.write("amr.max_level = 0\n")
        f.write("amr.blocking_factor = 8\n")
        f.write("amr.max_grid_size = 8\n")
        f.write("amr.plot_int = 1\n")
        f.write("geometry.coord_sys   = 0\n")
        f.write("geometry.is_periodic = 1 1 1\n")
        f.write("geometry.prob_lo = -0.5e-6 -0.5e-6 -0.5e-6\n")
        f.write("geometry.prob_hi = 0.5e-6 0.5e-6 0.5e-6\n")
        f.write("warpx.do_pml = 0\n")
        f.write("algo.charge_deposition = standard\n")
        f.write("algo.field_gathering = energy-conserving\n")
        f.write("warpx.cfl = 1.0\n")

        f.write("\nparticles.nspecies = {}\n".format(len(spec_names)))
        f.write("particles.species_names = {}\n".format(' '.join(spec_names)))
        f.write("particles.photon_species = {}\n".format(' '.join(spec_names)))

        f.write("\namr.plot_int = 50\n\n")

        for name in spec_names:
            f.write("{}.plot_species = 1\n".format(name))
            f.write("{}.plot_vars  = ux uy uz\n".format(name))

        f.write("\n")

        data = zip(spec_names, gamma_beta_list)
        for case in data:
            name = case[0]
            velx, vely ,velz = case[1]
            f.write("{}.species_type = photon\n".format(name))
            f.write('{}.injection_style = "SingleParticle"\n'.format(name))
            f.write("{}.single_particle_pos = {} {} {}\n".
                format(name, init_pos[0], init_pos[1], init_pos[2]))
            f.write("{}.single_particle_vel = {} {} {}\n".
                format(name, velx, vely, velz))
            f.write("{}.single_particle_weight = 1.0\n".format(name))
            f.write("\n".format(name))

def main():
    if (len(sys.argv) < 2):
        generate()
    else:
        check()

if __name__ == "__main__":
    main()
