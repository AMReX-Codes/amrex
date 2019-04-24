#! /usr/bin/env python

# This is a script that analyses the simulation results from
# the script `inputs.multi.rt`. This simulates a 3D periodic plasma wave.
# The electric field in the simulation is given (in theory) by:
# $$ E_x = \epsilon \,\frac{m_e c^2 k_x}{q_e}\sin(k_x x)\cos(k_y y)\cos(k_z z)\sin( \omega_p t)$$
# $$ E_y = \epsilon \,\frac{m_e c^2 k_y}{q_e}\cos(k_x x)\sin(k_y y)\cos(k_z z)\sin( \omega_p t)$$
# $$ E_z = \epsilon \,\frac{m_e c^2 k_z}{q_e}\cos(k_x x)\cos(k_y y)\sin(k_z z)\sin( \omega_p t)$$
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import yt
yt.funcs.mylog.setLevel(50)
import numpy as np
from scipy.constants import e, m_e, epsilon_0, c

# this will be the name of the plot file
fn = sys.argv[1]

# Parameters (these parameters must match the parameters in `inputs.multi.rt`)
epsilon = 0.01
n = 4.e24
n_osc_x = 2
n_osc_y = 2
n_osc_z = 2
xmin = -20e-6; xmax = 20.e-6; Nx = 64
ymin = -20e-6; ymax = 20.e-6; Ny = 64
zmin = -20e-6; zmax = 20.e-6; Nz = 64

# Wave vector of the wave
kx = 2.*np.pi*n_osc_x/(xmax-xmin)
ky = 2.*np.pi*n_osc_y/(ymax-ymin)
kz = 2.*np.pi*n_osc_z/(zmax-zmin)
# Plasma frequency
wp = np.sqrt((n*e**2)/(m_e*epsilon_0))

k = {'Ex':kx, 'Ey':ky, 'Ez':kz}
cos = {'Ex': (0,1,1), 'Ey':(1,0,1), 'Ez':(1,1,0)}

def get_contribution( is_cos, k ):
    du = (xmax-xmin)/Nx
    u = xmin + du*( 0.5 + np.arange(Nx) )
    if is_cos == 1:
        return( np.cos(k*u) )
    else:
        return( np.sin(k*u) )

def get_theoretical_field( field, t ):
    amplitude = epsilon * (m_e*c**2*k[field])/e * np.sin(wp*t)
    cos_flag = cos[field]
    x_contribution = get_contribution( cos_flag[0], kx )
    y_contribution = get_contribution( cos_flag[1], ky )
    z_contribution = get_contribution( cos_flag[2], kz )

    E = amplitude * x_contribution[:, np.newaxis, np.newaxis] \
                  * y_contribution[np.newaxis, :, np.newaxis] \
                  * z_contribution[np.newaxis, np.newaxis, :]

    return( E )

# Read the file
ds = yt.load(fn)

# Check that the particle selective output worked:
for species in ['electrons', 'positrons']:
    for field in ['particle_weight',
                  'particle_momentum_x',
                  'particle_Bz',
                  'particle_Ey']:
        assert (species, field) in ds.field_list
    for field in ['particle_momentum_y',
                  'particle_momentum_z',
                  'particle_Bx',
                  'particle_By',
                  'particle_Ex',
                  'particle_Ez']:
        assert (species, field) not in ds.field_list


t0 = ds.current_time.to_ndarray().mean()
data = ds.covering_grid(level=0, left_edge=ds.domain_left_edge,
                                    dims=ds.domain_dimensions)

# Check the validity of the fields
overall_max_error = 0
for field in ['Ex', 'Ey', 'Ez']:
    E_sim = data[field].to_ndarray()
    E_th = get_theoretical_field(field, t0)
    max_error = abs(E_sim-E_th).max()/abs(E_th).max()
    print('%s: Max error: %.2e' %(field,max_error))
    overall_max_error = max( overall_max_error, max_error )

# Plot the last field from the loop (Ez at iteration 40)
plt.subplot2grid( (1,2), (0,0) )
plt.imshow( E_sim[:,:,32] )
plt.colorbar()
plt.title('Ez, last iteration\n(simulation)')
plt.subplot2grid( (1,2), (0,1) )
plt.imshow( E_th[:,:,32] )
plt.colorbar()
plt.title('Ez, last iteration\n(theory)')
plt.tight_layout()
plt.savefig('langmuir_multi_analysis.png')

# Automatically check the validity
assert overall_max_error < 0.035
