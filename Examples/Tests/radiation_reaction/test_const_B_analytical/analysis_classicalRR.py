#! /usr/bin/env python

# Copyright 2019 Luca Fedeli, Maxence Thevenet, Remi Lehe
#
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL


# This script contains few simple tests for the radiation reaction pusher
# It initializes an electron or a positron with normalized momentum in different
# directions, propagating in a static magnetic field (along [2/7,3/7,6/7]).
# If the initial momentum is perpendicular to the the field, there is a very
# simple analytical formula for the gamma factor:
# gamma(t) =  coth(t/tc - C)
# where tc = 1./(t0 * omega_c * omega_c)
# and omega_c =  qB/m , t0 = (2/3) * (r_e / c)
# and r_e = classical electron radius
# C is chosen so that gamma(0) = initial_gamma
# If the initial momentum is parallel to the field, it should not change.
# This test checks these predictions with a tolerance of 5%.
# If the script is run without a command line argument, it regenerates a new
# inputfile according to the initial conditions listed below.
# For detailed references see:
# 1) M. Tamburini. PhD Thesis. University of Pisa (2011)
#   https://etd.adm.unipi.it/t/etd-11232011-111059/
# 2) H. Spohn Europhysics Letters 50 287 (2000)
#   https://arxiv.org/abs/physics/9911027
# 3) H. Spohn, Dynamics of charged particles and their radiation field
#   (Cambridge University Press, Cambridge, 2004)

import numpy as np
import sys
import yt

#Input filename
inputname = "inputs"
#________________________________________

#Physical constants
c = 299792458.
m_e = 9.1093837015e-31
q_0 = 1.602176634e-19
classical_electron_radius = 2.81794e-15
reference_length = 1.0e-6
very_small_dot_product = 1.0e-4
very_small_weight = 1.0e-8
#________________________________________

#Sim box data
sim_size = 0.8e-6
resolution = 64
steps = 64
init_pos = np.array([0.,0.,0.])
#________________________________________

#Momentum vals
p_aux_0 = np.array([2.,3.,6.])
p_aux_1 = np.array([1,0,0])
p_aux_2 = np.array([0,1,0])
Q, _ = np.linalg.qr(np.column_stack( [p_aux_0, p_aux_1, p_aux_2] )) #Gram-Schmidt
p_0 = -Q[:,0]
p_1 = -Q[:,1]
p_2 = -Q[:,2]

p_vals = [50,200,1000]
#________________________________________

#Field val
B_val_norm = 300
B_val = B_val_norm*m_e * 2.0 * np.pi * c /q_0/reference_length
B = p_0 * B_val
#________________________________________

#Tolerance
tol = 0.05
#________________________________________

#tau_c
omega_c = q_0*B_val/m_e
t0 = (2./3.)*classical_electron_radius/c
tau_c = 1.0/omega_c/omega_c/t0
#________________________________________


#Simulation case struct
class sim_case:
    def __init__(self, _name, _init_mom, _type):
        self.name = _name
        self.init_mom = _init_mom
        self.type = _type
#________________________________________

#All cases
cases = [
    sim_case("ele_para0", p_0 * p_vals[2], "-q_e"),
    sim_case("ele_perp0", p_1 * p_vals[0], "-q_e"),
    sim_case("ele_perp1", p_2 * p_vals[1], "-q_e"),
    sim_case("ele_perp2", p_1 * p_vals[2], "-q_e"),
    sim_case("pos_perp2", p_1 * p_vals[2], "q_e"),
]
#________________________________________

#Auxiliary functions
def gamma(p) :
    return np.sqrt(1.0 + np.dot(p,p))

def exp_res(cc, time):
    if np.all(np.linalg.norm(np.cross(cc.init_mom, B)) < very_small_dot_product):
        return gamma(cc.init_mom)
    else :
        tt = time/tau_c
        g0 = gamma(cc.init_mom)
        C = -0.5 * np.log((g0+1)/(g0-1))
        return 1.0/np.tanh(tt - C)
#________________________________________


def check():
    filename = sys.argv[1]
    data_set_end = yt.load(filename)

    sim_time = data_set_end.current_time.to_value()

    #simulation results
    all_data =  data_set_end.all_data()
    spec_names = [cc.name for cc in cases]

    #All momenta
    res_mom = np.array([np.array([
        all_data[sp, 'particle_momentum_x'].v[0],
        all_data[sp, 'particle_momentum_y'].v[0],
        all_data[sp, 'particle_momentum_z'].v[0]])
        for sp in spec_names])

    for cc in zip(cases, res_mom):
        init_gamma = gamma(cc[0].init_mom)
        end_gamma = gamma(cc[1]/m_e/c)
        exp_gamma = exp_res(cc[0], sim_time)
        assert(np.abs(end_gamma-exp_gamma)/exp_gamma < tol)

def generate():

    with open(inputname,'w') as f:
        f.write("#Automatically generated inputfile\n")
        f.write("#Run check.py without arguments to regenerate\n")
        f.write("#\n\n")
        f.write("max_step = {}\n".format(steps))
        f.write("amr.n_cell = {} {} {}\n".format(resolution, resolution, resolution))
        f.write("amr.max_level = 0\n")
        f.write("amr.blocking_factor = 32\n")
        f.write("amr.max_grid_size = 64\n")
        f.write("geometry.coord_sys   = 0\n")
        f.write("geometry.is_periodic = 1 1 1\n")
        f.write("geometry.prob_lo = {} {} {}\n".format(-sim_size, -sim_size, -sim_size))
        f.write("geometry.prob_hi = {} {} {}\n".format(sim_size, sim_size, sim_size))
        f.write("warpx.do_pml = 0\n")
        f.write("algo.charge_deposition = standard\n")
        f.write("algo.field_gathering = energy-conserving\n")
        f.write("warpx.cfl = 1.0\n")
        f.write("warpx.serialize_ics = 1\n")

        f.write("\nparticles.nspecies = {}\n".format(len(cases)))
        f.write("particles.species_names = ")
        for cc in cases:
            f.write(" {}".format(cc.name))
        f.write("\n")

        f.write("\namr.plot_int = {}\n\n".format(steps))

        f.write("warpx.fields_to_plot = rho\n\n")

        for cc in cases:
            f.write("{}.charge = {}\n".format(cc.name, cc.type))
            f.write("{}.mass = m_e\n".format(cc.name))
            f.write('{}.injection_style = "SingleParticle"\n'.format(cc.name))
            f.write("{}.single_particle_pos = {} {} {}\n".
                format(cc.name, init_pos[0], init_pos[1], init_pos[2]))
            f.write("{}.single_particle_vel = {} {} {}\n".
                format(cc.name, cc.init_mom[0], cc.init_mom[1], cc.init_mom[2]))
            f.write("{}.single_particle_weight = {}\n".format(cc.name, very_small_weight))
            f.write("{}.do_classical_radiation_reaction = 1\n".format(cc.name))
            f.write("\n")

        f.write("warpx.B_external_particle = {} {} {}\n".format(*B))

def main():
    if (len(sys.argv) < 2):
        generate()
    else:
        check()

if __name__ == "__main__":
    main()
