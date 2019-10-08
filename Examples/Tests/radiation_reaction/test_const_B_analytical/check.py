#! /usr/bin/env python3

# This script contains few simple tests for the radiation reaction pusher
# It initializes an electron or a positron with normalized momentum in different
# directions, propagating in a static magnetic field (along z).
# If the initial momentum is perpendicular to the the field, there is a very
# simple analytical formula for the gamma factor:
# gamma(t) =  coth(t/tc - C)
# where tc = 1./(t0 * omega_c * omega_c)
# and omega_c =  qB/m , t0 = (2/3) * (r_e / c)
# and r_e = classical electron radius
# C is chosen so that gamma(0) = initial_gamma
# If the initial momentum is parallel to the field, it should not change.
# This module tests these predictions with a tolerance of 5%
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
import itertools as itr


#Input filename
inputname = "inputs"
#________________________________________

#Physical constants
c = 299792458.
m_e = 9.1093837015e-31
q_0 = 1.602176634e-19
classical_electron_radius = 2.81794e-15
reference_length = 1.0e-6
small = 1.0e-100
#________________________________________

#Sim box data
sim_size = 1.0e-6
resolution = 128
steps = 128
init_pos = np.array([0.,0.,0.])
#________________________________________

#Momentum vals
p_vals = [10, 20, 50, 100, 200, 500, 1000]
#________________________________________

#Field val
B_val_norm = 300
B_val = B_val_norm*m_e * 2.0 * np.pi * c /q_0/reference_length
B = np.array([0, 0, B_val])
#________________________________________

#Dirs
dirs = (np.array([1.,0.,0.]), np.array([0.,1.,0.]), np.array([0.,0.,1.]))
#________________________________________

#Types
types = ["q_e", "-q_e"]
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
    sim_case("s" + str(ii),
        vv[0]*vv[1], # momentum
        vv[2]) #type
        for ii, vv in enumerate(list(itr.product(p_vals, dirs, types)))
]
#________________________________________

#Auxiliary functions
def gamma(p) :
    return np.sqrt(1.0 + np.dot(p,p))

def exp_res(cc, time):
    if np.all(np.cross(cc.init_mom, B) == 0.):
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
        assert(np.all(np.abs(end_gamma-exp_gamma)/exp_gamma < tol))

def generate():

    with open(inputname,'w') as f:
        f.write("#Automatically generated inputfile\n")
        f.write("#Run check.py without arguments to regenerate\n")
        f.write("#\n\n")
        f.write("max_step = {}\n".format(steps))
        f.write("amr.n_cell = {} {} {}\n".format(resolution, resolution, resolution))
        f.write("amr.max_level = 0\n")
        f.write("amr.blocking_factor = 8\n")
        f.write("amr.max_grid_size = 8\n")
        f.write("geometry.coord_sys   = 0\n")
        f.write("geometry.is_periodic = 1 1 1\n")
        f.write("geometry.prob_lo = {} {} {}\n".format(-sim_size, -sim_size, -sim_size))
        f.write("geometry.prob_hi = {} {} {}\n".format(sim_size, sim_size, sim_size))
        f.write("warpx.do_pml = 0\n")
        f.write("algo.charge_deposition = standard\n")
        f.write("algo.field_gathering = standard\n")
        f.write("warpx.cfl = 1.0\n")

        f.write("\nparticles.nspecies = {}\n".format(len(cases)))
        f.write("particles.species_names = ")
        for cc in cases:
            f.write(" {}".format(cc.name))
        f.write("\n")

        f.write("\namr.plot_int = {}\n\n".format(steps))

        f.write("warpx.plot_raw_fields = 0\n\n")

        for cc in cases:
            f.write("{}.charge = {}\n".format(cc.name, cc.type))
            f.write("{}.mass = m_e\n".format(cc.name))
            f.write('{}.injection_style = "SingleParticle"\n'.format(cc.name))
            f.write("{}.single_particle_pos = {} {} {}\n".
                format(cc.name, init_pos[0], init_pos[1], init_pos[2]))
            f.write("{}.single_particle_vel = {} {} {}\n".
                format(cc.name, cc.init_mom[0], cc.init_mom[1], cc.init_mom[2]))
            f.write("{}.single_particle_weight = {}\n".format(cc.name, small))
            f.write("{}.do_classical_radiation_reaction = 1 \n".format(cc.name))
            f.write("\n")

        f.write("warpx.B_external = {} {} {}\n".format(*B))


def main():
    if (len(sys.argv) < 2):
        generate()
    else:
        check()

if __name__ == "__main__":
    main()
