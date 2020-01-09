# This is a script that analyses the multimode simulation results.
# This simulates a RZ multimode periodic plasma wave.
# The electric field from the simulation is compared to the analytic value

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pywarpx import picmi

constants = picmi.constants

nr = 64
nz = 200

rmin =  0.e0
zmin =  0.e0
rmax = +20.e-6
zmax = +40.e-6

# Parameters describing particle distribution
density = 2.e24
epsilon0 = 0.001*constants.c
epsilon1 = 0.001*constants.c
epsilon2 = 0.001*constants.c
w0 = 5.e-6
n_osc_z = 3

# Wave vector of the wave
k0 = 2.*np.pi*n_osc_z/(zmax - zmin)

# Plasma frequency
wp = np.sqrt((density*constants.q_e**2)/(constants.m_e*constants.ep0))
kp = wp/constants.c

uniform_plasma = picmi.UniformDistribution(density = density,
                                           upper_bound = [+18e-6, +18e-6, None],
                                           directed_velocity = [0., 0., 0.])

momentum_expressions = ["""+ epsilon0/kp*2*x/w0**2*exp(-(x**2+y**2)/w0**2)*sin(k0*z)
                           - epsilon1/kp*2/w0*exp(-(x**2+y**2)/w0**2)*sin(k0*z)
                           + epsilon1/kp*4*x**2/w0**3*exp(-(x**2+y**2)/w0**2)*sin(k0*z)
                           - epsilon2/kp*8*x/w0**2*exp(-(x**2+y**2)/w0**2)*sin(k0*z)
                           + epsilon2/kp*8*x*(x**2-y**2)/w0**4*exp(-(x**2+y**2)/w0**2)*sin(k0*z)""",
                        """+ epsilon0/kp*2*y/w0**2*exp(-(x**2+y**2)/w0**2)*sin(k0*z)
                           + epsilon1/kp*4*x*y/w0**3*exp(-(x**2+y**2)/w0**2)*sin(k0*z)
                           + epsilon2/kp*8*y/w0**2*exp(-(x**2+y**2)/w0**2)*sin(k0*z)
                           + epsilon2/kp*8*y*(x**2-y**2)/w0**4*exp(-(x**2+y**2)/w0**2)*sin(k0*z)""",
                        """- epsilon0/kp*k0*exp(-(x**2+y**2)/w0**2)*cos(k0*z)
                           - epsilon1/kp*k0*2*x/w0*exp(-(x**2+y**2)/w0**2)*cos(k0*z)
                           - epsilon2/kp*k0*4*(x**2-y**2)/w0**2*exp(-(x**2+y**2)/w0**2)*cos(k0*z)"""]

analytic_plasma = picmi.AnalyticDistribution(density_expression = density,
                                             upper_bound = [+18e-6, +18e-6, None],
                                             epsilon0 = epsilon0,
                                             epsilon1 = epsilon1,
                                             epsilon2 = epsilon2,
                                             kp = kp,
                                             k0 = k0,
                                             w0 = w0,
                                             momentum_expressions = momentum_expressions)

electrons = picmi.Species(particle_type='electron', name='electrons', initial_distribution=analytic_plasma)
protons = picmi.Species(particle_type='proton', name='protons', initial_distribution=uniform_plasma)

grid = picmi.CylindricalGrid(number_of_cells = [nr, nz],
                             n_azimuthal_modes = 3,
                             lower_bound = [rmin, zmin],
                             upper_bound = [rmax, zmax],
                             lower_boundary_conditions = ['dirichlet', 'periodic'],
                             upper_boundary_conditions = ['dirichlet', 'periodic'],
                             moving_window_zvelocity = 0.,
                             warpx_max_grid_size=64)

solver = picmi.ElectromagneticSolver(grid=grid, cfl=1.)

sim = picmi.Simulation(solver = solver,
                       max_steps = 40,
                       verbose = 1,
                       warpx_plot_int = 40,
                       warpx_current_deposition_algo = 'esirkepov',
                       warpx_field_gathering_algo = 'energy-conserving',
                       warpx_particle_pusher_algo = 'boris')

sim.add_species(electrons, layout=picmi.GriddedLayout(n_macroparticle_per_cell=[2,16,2], grid=grid))
sim.add_species(protons, layout=picmi.GriddedLayout(n_macroparticle_per_cell=[2,16,2], grid=grid))

# write_inputs will create an inputs file that can be used to run
# with the compiled version.
#sim.write_input_file(file_name='inputsrz_from_PICMI')

# Alternatively, sim.step will run WarpX, controlling it from Python
sim.step()


# Below is WarpX specific code to check the results.
import pywarpx
from pywarpx.fields import *

def calcEr( z, r, k0, w0, wp, t, epsilons) :
    """
    Return the radial electric field as an array
    of the same length as z and r, in the half-plane theta=0
    """
    Er_array = (
        epsilons[0] * constants.m_e*constants.c/constants.q_e * 2*r/w0**2 *
            np.exp( -r**2/w0**2 ) * np.sin( k0*z ) * np.sin( wp*t )
        - epsilons[1] * constants.m_e*constants.c/constants.q_e * 2/w0 *
            np.exp( -r**2/w0**2 ) * np.sin( k0*z ) * np.sin( wp*t )
        + epsilons[1] * constants.m_e*constants.c/constants.q_e * 4*r**2/w0**3 *
            np.exp( -r**2/w0**2 ) * np.sin( k0*z ) * np.sin( wp*t )
        - epsilons[2] * constants.m_e*constants.c/constants.q_e * 8*r/w0**2 *
            np.exp( -r**2/w0**2 ) * np.sin( k0*z ) * np.sin( wp*t )
        + epsilons[2] * constants.m_e*constants.c/constants.q_e * 8*r**3/w0**4 *
            np.exp( -r**2/w0**2 ) * np.sin( k0*z ) * np.sin( wp*t ))
    return( Er_array )

def calcEz( z, r, k0, w0, wp, t, epsilons) :
    """
    Return the longitudinal electric field as an array
    of the same length as z and r, in the half-plane theta=0
    """
    Ez_array = (
        - epsilons[0] * constants.m_e*constants.c/constants.q_e * k0 *
            np.exp( -r**2/w0**2 ) * np.cos( k0*z ) * np.sin( wp*t )
        - epsilons[1] * constants.m_e*constants.c/constants.q_e * k0 * 2*r/w0 *
            np.exp( -r**2/w0**2 ) * np.cos( k0*z ) * np.sin( wp*t )
        - epsilons[2] * constants.m_e*constants.c/constants.q_e * k0 * 4*r**2/w0**2 *
            np.exp( -r**2/w0**2 ) * np.cos( k0*z ) * np.sin( wp*t ))
    return( Ez_array )

# Get node centered coordinates
dr = (rmax - rmin)/nr
dz = (zmax - zmin)/nz
coords = np.indices([nr+1, nz+1], 'd')
rr = rmin + coords[0]*dr
zz = zmin + coords[1]*dz

# Current time of the simulation
t0 = pywarpx._libwarpx.libwarpx.warpx_gett_new(0)

# Get the raw field data. Note that these are the real and imaginary
# parts of the fields for each azimuthal mode.
Ex_sim_modes = ExWrapper()[...]
Ez_sim_modes = EzWrapper()[...]

# Sum the real components to get the field along x-axis (theta = 0)
Er_sim = Ex_sim_modes[:,:,0] + np.sum(Ex_sim_modes[:,:,1::2], axis=2)
Ez_sim = Ez_sim_modes[:,:,0] + np.sum(Ez_sim_modes[:,:,1::2], axis=2)

# The analytical solutions
Er_th = calcEr(zz[:-1,:], rr[:-1,:] + dr/2., k0, w0, wp, t0, [epsilon0, epsilon1, epsilon2])
Ez_th = calcEz(zz[:,:-1] + dz/2., rr[:,:-1], k0, w0, wp, t0, [epsilon0, epsilon1, epsilon2])

max_error_Er = abs(Er_sim - Er_th).max()/abs(Er_th).max()
max_error_Ez = abs(Ez_sim - Ez_th).max()/abs(Ez_th).max()
print("Max error Er %e"%max_error_Er)
print("Max error Ez %e"%max_error_Ez)

# Plot the last field from the loop (Ez at iteration 40)
plt.subplot2grid( (1,2), (0,0) )
plt.imshow( Ez_sim )
plt.colorbar()
plt.title('Ez, last iteration\n(simulation)')
plt.subplot2grid( (1,2), (0,1) )
plt.imshow( Ez_th )
plt.colorbar()
plt.title('Ez, last iteration\n(theory)')
plt.tight_layout()
plt.savefig('langmuir_multi_rz_multimode_analysis.png')

assert max(max_error_Er, max_error_Ez) < 0.02
