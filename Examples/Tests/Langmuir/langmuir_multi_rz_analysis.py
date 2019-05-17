#! /usr/bin/env python

# This is a script that analyses the simulation results from
# the script `inputs.multi.rz.rt`. This simulates a RZ periodic plasma wave.
# The electric field in the simulation is given (in theory) by:
# $$ E_r = -\partial_r \phi = \epsilon \,\frac{mc^2}{e}\frac{2\,r}{w_0^2} \exp\left(-\frac{r^2}{w_0^2}\right) \sin(k_0 z) \sin(\omega_p t)
# $$ E_z = -\partial_z \phi = - \epsilon \,\frac{mc^2}{e} k_0 \exp\left(-\frac{r^2}{w_0^2}\right) \cos(k_0 z) \sin(\omega_p t)
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

# Parameters (these parameters must match the parameters in `inputs.multi.rz.rt`)
epsilon = 0.01
n = 2.e24
w0 = 5.e-6
n_osc_z = 2
rmin =   0e-6; rmax = 20.e-6; Nr = 64
zmin = -20e-6; zmax = 20.e-6; Nz = 128

# Wave vector of the wave
k0 = 2.*np.pi*n_osc_z/(zmax-zmin)
# Plasma frequency
wp = np.sqrt((n*e**2)/(m_e*epsilon_0))
kp = wp/c

def Er( z, r, epsilon, k0, w0, wp, t) :
    """
    Return the radial electric field as an array
    of the same length as z and r, in the half-plane theta=0
    """
    Er_array = \
        epsilon * m_e*c**2/e * 2*r/w0**2 * \
            np.exp( -r**2/w0**2 ) * np.sin( k0*z ) * np.sin( wp*t )
    return( Er_array )

def Ez( z, r, epsilon, k0, w0, wp, t) :
    """
    Return the longitudinal electric field as an array
    of the same length as z and r, in the half-plane theta=0
    """
    Ez_array = \
        - epsilon * m_e*c**2/e * k0 * \
            np.exp( -r**2/w0**2 ) * np.cos( k0*z ) * np.sin( wp*t )
    return( Ez_array )

# Read the file
ds = yt.load(fn)
t0 = ds.current_time.to_ndarray().mean()
data = ds.covering_grid(level=0, left_edge=ds.domain_left_edge,
                        dims=ds.domain_dimensions)

# Get cell centered coordinates
dr = (rmax - rmin)/Nr
dz = (zmax - zmin)/Nz
coords = np.indices([Nr, Nz],'d')
rr = rmin + (coords[0] + 0.5)*dr
zz = zmin + (coords[1] + 0.5)*dz

# Check the validity of the fields
overall_max_error = 0
Er_sim = data['Ex'].to_ndarray()[:,:,0]
Er_th = Er(zz, rr, epsilon, k0, w0, wp, t0)
max_error = abs(Er_sim-Er_th).max()/abs(Er_th).max()
print('Er: Max error: %.2e' %(max_error))
overall_max_error = max( overall_max_error, max_error )

Ez_sim = data['Ez'].to_ndarray()[:,:,0]
Ez_th = Ez(zz, rr, epsilon, k0, w0, wp, t0)
max_error = abs(Ez_sim-Ez_th).max()/abs(Ez_th).max()
print('Ez: Max error: %.2e' %(max_error))
overall_max_error = max( overall_max_error, max_error )

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
plt.savefig('langmuir_multi_rz_analysis.png')

# Automatically check the validity
assert overall_max_error < 0.04
