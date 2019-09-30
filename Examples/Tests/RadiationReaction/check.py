#! /usr/bin/env python3

# This script checks if the Radiation Reaction pusher is working properly.
# It initializes particles with different momenta in a static electromagnetic
# field and it checks if the final momenta and positions are close enough to
# the expected ones. The expected final positions and momenta are calculated
# using a solver (validated separately)

# If the script is run without a command line argument, it regenerates a new
# inputfile according to the initial conditions listed below.

import numpy as np
import sys
sys.path.append('./aux')
import solver
import test_solver

import yt


#Physical constants
c = 299792458.
m_e = 9.1093837015e-31
q_0 = 1.602176634e-19
epsi_0 =  8.8541878128e-12
mu_0 = 4.0*np.pi*1e-7
reference_length = 1.0e-6
#________________________________________

#Conversion factors
norm_to_si = {
"E" : (m_e * 2.0 * np.pi * c /q_0/reference_length)/np.sqrt(4.0*np.pi*epsi_0),
"B" : (m_e * 2.0 * np.pi * c /q_0/reference_length)/np.sqrt(4.0*np.pi/mu_0)
}
si_to_norm = {
"x" : 1.0/(2.0*np.pi/reference_length),
"t" : 1.0/(2.0*np.pi*c/reference_length),
"p" : 1.0/m_e/c
}
#________________________________________


#Initial conditions
EVAL = np.array([4321, -1543, 3151])
BVAL = np.array([-2971, 3967, -4332])
EVAL_SI =  norm_to_si["E"] * EVAL
BVAL_SI = norm_to_si["B"] * BVAL
E_NORM = lambda pos, t : EVAL
B_NORM = lambda pos, t : BVAL
E_SI = lambda pos, t : EVAL_SI
B_SI = lambda pos, t : BVAL_SI
init_mom = np.array([34, -12, 81])
init_pos= np.array([0.0, 0.0, 0.0])
#________________________________________

#Input filename
inputname = "inputs"
#________________________________________

# This function reads the WarpX plotfile given as the first command-line
# argument, and check if the position and momenta of particles agrees with theory.
def check():
    filename = sys.argv[1]
    data_set_end = yt.load(filename)

    sim_time = data_set_end.current_time.to_value()
    print(sim_time)

    #simulation results
    all_data =  data_set_end.all_data()
    spec_names = ("ele", "pos")
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

    sol_ele = solver.solve(
        init_pos*si_to_norm["x"],
        init_mom,
        0.0,
        sim_time*si_to_norm["t"],
        20000, E_NORM, B_NORM, -1, reference_length)

    sol_pos = solver.solve(
        init_pos*si_to_norm["x"],
        init_mom,
        0.0,
        sim_time*si_to_norm["t"],
        20000, E_NORM, B_NORM, 1, reference_length)

    print(sol_ele[-1,0], sol_ele[-1,1], sol_ele[-1,2], sol_pos[-1,0], sol_pos[-1,1], sol_pos[-1,2])
    print(sol_ele[-1,3], sol_ele[-1,4], sol_ele[-1,5], sol_pos[-1,3], sol_pos[-1,4], sol_pos[-1,5])
    print(np.array(res_pos)*si_to_norm["x"])
    print(np.array(res_mom)*si_to_norm["p"])

# This function generates the input file to test the radiation reaction pusher
def generate():
    # cheks if the solver is working properly
    test_solver.do_test()

    with open(inputname,'w') as f:
        f.write("#Automatically generated inputfile\n")
        f.write("#Run check.py without arguments to regenerate\n")
        f.write("#\n\n")
        f.write("max_step = 100\n")
        f.write("amr.n_cell = 128 128 128\n")
        f.write("amr.max_level = 0\n")
        f.write("amr.blocking_factor = 8\n")
        f.write("amr.max_grid_size = 8\n")
        f.write("amr.plot_int = 1\n")
        f.write("geometry.coord_sys   = 0\n")
        f.write("geometry.is_periodic = 1 1 1\n")
        f.write("geometry.prob_lo = -1e-6 -1e-6 -1e-6\n")
        f.write("geometry.prob_hi = 1e-6 1e-6 1e-6\n")
        f.write("warpx.do_pml = 0\n")
        f.write("algo.charge_deposition = standard\n")
        f.write("algo.field_gathering = standard\n")
        f.write("warpx.cfl = 1.0\n")

        f.write("\nparticles.nspecies = 2\n")
        f.write("particles.species_names = ele pos")

        f.write("\namr.plot_int = 200\n\n")

        for name in ("ele", "pos"):
            f.write("{}.plot_species = 1\n".format(name))
            f.write("{}.plot_vars  = ux uy uz\n".format(name))

        f.write("\n")

        for case in  zip(("ele", "pos"), ("-q_e", "q_e")):
            name = case[0]
            velx, vely ,velz = init_mom
            f.write("{}.charge = {}\n".format(name, case[1]))
            f.write("{}.mass = m_e\n".format(name))
            f.write('{}.injection_style = "SingleParticle"\n'.format(name))
            f.write("{}.single_particle_pos = {} {} {}\n".
                format(name, init_pos[0], init_pos[1], init_pos[2]))
            f.write("{}.single_particle_vel = {} {} {}\n".
                format(name, velx, vely, velz))

            # An extremely small weigth (ideally I would like a test particle)
            f.write("{}.single_particle_weight = 1.0e-64\n".format(name))
            f.write("{}.do_classical_radiation_reaction = 1 \n".format(name))
            f.write("\n")

        f.write("warpx.E_external = {} {} {}\n".format(*EVAL_SI))
        f.write("warpx.B_external = {} {} {}\n".format(*BVAL_SI))

def main():
    if (len(sys.argv) < 2):
        generate()
    else:
        check()

if __name__ == "__main__":
    main()
