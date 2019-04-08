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

# Parameters (these parameters must match the parameters in `inputs.multi.rt`)
epsilon = 0.001
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

epsilons = [epsilon, 0., 0.]

def Jx( z, r, epsilons, k0, w0, wp, t) :
    return -n*e*c*epsilons[0]/kp*2*r/w0**2*np.exp(-(r**2)/w0**2)*np.sin(k0*z)*np.cos(wp*t)

def Jz( z, r, epsilons, k0, w0, wp, t) :
    return +n*e*c*epsilons[0]/kp*k0*np.exp(-(r**2)/w0**2)*np.cos(k0*z)*np.cos(wp*t)

def Er( z, r, epsilons, k0, w0, wp, t) :
    """
    Return the radial electric field as an array
    of the same length as z and r, in the half-plane theta=0
    """
    Er_array = \
        epsilons[0] * m_e*c**2/e * 2*r/w0**2 * \
            np.exp( -r**2/w0**2 ) * np.sin( k0*z ) * np.sin( wp*t ) \
        - epsilons[1] * m_e*c**2/e * 2/w0 * \
            np.exp( -r**2/w0**2 ) * np.sin( k0*z ) * np.sin( wp*t ) \
        + epsilons[1] * m_e*c**2/e * 4*r**2/w0**3 * \
            np.exp( -r**2/w0**2 ) * np.sin( k0*z ) * np.sin( wp*t ) \
        - epsilons[2] * m_e*c**2/e * 8*r/w0**2 * \
            np.exp( -r**2/w0**2 ) * np.sin( k0*z ) * np.sin( wp*t ) \
        + epsilons[2] * m_e*c**2/e * 8*r**3/w0**4 * \
            np.exp( -r**2/w0**2 ) * np.sin( k0*z ) * np.sin( wp*t )
    return( Er_array )

def Ez( z, r, epsilons, k0, w0, wp, t) :
    """
    Return the longitudinal electric field as an array
    of the same length as z and r, in the half-plane theta=0
    """
    Ez_array = \
        - epsilons[0] * m_e*c**2/e * k0 * \
            np.exp( -r**2/w0**2 ) * np.cos( k0*z ) * np.sin( wp*t ) \
        - epsilons[1] * m_e*c**2/e * k0 * 2*r/w0 * \
            np.exp( -r**2/w0**2 ) * np.cos( k0*z ) * np.sin( wp*t ) \
        - epsilons[2] * m_e*c**2/e * k0 * 4*r**2/w0**2 * \
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
Er_th = Er(zz, rr, epsilons, k0, w0, wp, t0)
max_error = abs(Er_sim-Er_th).max()/abs(Er_th).max()
print('Er: Max error: %.2e' %(max_error))
overall_max_error = max( overall_max_error, max_error )

Ez_sim = data['Ez'].to_ndarray()[:,:,0]
Ez_th = Ez(zz, rr, epsilons, k0, w0, wp, t0)
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


"""
# Extra plots
import matplotlib.cm as cm

# Side by side plot of Er
fig, ax = plt.subplots(1,2)
im = ax[0].imshow(Er_sim, interpolation='bilinear', origin='lower', cmap=cm.RdBu, extent=(-20, 20, 0, 20))
fig.colorbar(im, ax=ax[0], orientation='vertical', shrink=0.6)
ax[0].set_title('Er sim')
im = ax[1].imshow(Er_th, interpolation='bilinear', origin='lower', cmap=cm.RdBu, extent=(-20, 20, 0, 20))
fig.colorbar(im, ax=ax[1], orientation='vertical', shrink=0.6)
ax[1].set_title('Er th')
fig.show()

# Side by side plot of Ez
fig, ax = plt.subplots(1,2)
im = ax[0].imshow(Ez_sim, interpolation='bilinear', origin='lower', cmap=cm.RdBu, extent=(-20, 20, 0, 20))
fig.colorbar(im, ax=ax[0], orientation='vertical', shrink=0.6)
ax[0].set_title('Ez sim')
im = ax[1].imshow(Ez_th, interpolation='bilinear', origin='lower', cmap=cm.RdBu, extent=(-20, 20, 0, 20))
fig.colorbar(im, ax=ax[1], orientation='vertical', shrink=0.6)
ax[1].set_title('Ez th')
fig.show()

# Differences of Er and Ez
fig, ax = plt.subplots(1,2)
im = ax[0].imshow((Er_sim - Er_th)/Er_sim.max(), interpolation='bilinear', origin='lower', cmap=cm.RdBu, extent=(-20, 20, 0, 20))
fig.colorbar(im, ax=ax[0], orientation='vertical', shrink=0.5)
ax[0].set_title('Er difference')
im = ax[1].imshow((Ez_sim - Ez_th)/Ez_sim.max(), interpolation='bilinear', origin='lower', cmap=cm.RdBu, extent=(-20, 20, 0, 20))
fig.colorbar(im, ax=ax[1], orientation='vertical', shrink=0.5)
ax[1].set_title('Ez difference')
fig.show()

# Er along z
fig, ax = plt.subplots()
ax.plot(zz[0,:]*1.e6, Er_sim[2,:], 'r')
ax.plot(zz[0,:]*1.e6, Er_th[2,:], 'b')
fig.show()

# Er along radius
fig, ax = plt.subplots()
ax.plot(rr[:,48]*1.e6, Er_sim[:,48], 'r')
ax.plot(rr[:,48]*1.e6, Er_th[:,48], 'b')
fig.show()

# Er difference along radius
fig, ax = plt.subplots()
ax.plot(rr[:,48]*1.e6, (Er_sim[:,48] - Er_th[:,48])/Er_sim.max(), 'r')
fig.show()

# Ez along radius
fig, ax = plt.subplots()
ax.plot(rr[:,48]*1.e6, Ez_sim[:,48], 'r')
ax.plot(rr[:,48]*1.e6, Ez_th[:,48], 'b')
fig.show()

# Ez difference along radius
fig, ax = plt.subplots()
ax.plot(rr[:,48]*1.e6, (Ez_sim[:,48] - Ez_th[:,48])/Ez_sim.max(), 'r')
fig.show()

# Ez along z
fig, ax = plt.subplots()
ax.plot(zz[0,:]*1.e6, Ez_sim[2,:], 'r')
ax.plot(zz[0,:]*1.e6, Ez_th[2,:], 'b')
fig.show()

# Fetch Jr
Jr_sim = data['Jx'].to_ndarray()[:,:,0]
Jr_th = Jx(zz, rr, epsilons, k0, w0, wp, t0)
max_error = abs(Jr_sim-Jr_th).max()/abs(Jr_th).max()
print('Jr: Max error: %.2e' %(max_error))

# Fetch Jz
Jz_sim = data['Jz'].to_ndarray()[:,:,0]
Jz_th = Jz(zz, rr, epsilons, k0, w0, wp, t0)
max_error = abs(Jz_sim-Jz_th).max()/abs(Jz_th).max()
print('Jz: Max error: %.2e' %(max_error))

# Plot Jr along radius
fig, ax = plt.subplots()
ax.plot(rr[:,48]*1.e6, Jr_sim[:,48], 'r')
ax.plot(rr[:,48]*1.e6, Jr_th[:,48], 'b')
fig.show()

"""
