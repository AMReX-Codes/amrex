# This is a script that analyses the simulation results from
# the script `inputs.multi.rt`. This simulates a 3D periodic plasma wave.
# The electric field in the simulation is given (in theory) by:
# $$ E_x = -\epsilon \,\frac{m_e c^2 k_x}{q_e}\cos(k_x x)\sin(k_y y)\sin(k_z z)\sin( \omega_p t)$$
# $$ E_y = -\epsilon \,\frac{m_e c^2 k_y}{q_e}\sin(k_x x)\cos(k_y y)\sin(k_z z)\sin( \omega_p t)$$
# $$ E_z = -\epsilon \,\frac{m_e c^2 k_z}{q_e}\sin(k_x x)\sin(k_y y)\cos(k_z z)\sin( \omega_p t)$$

import yt
import numpy as np
from scipy.constants import e, m_e, epsilon_0, c
import matplotlib.pyplot as plt
from tqdm import tqdm

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
sin = {'Ex': (0,1,1), 'Ey':(1,0,1), 'Ez':(1,1,0)}

def get_contribution( is_sin, k ):
    du = (xmax-xmin)/Nx
    u = xmin + du*( 0.5 + np.arange(Nx) )
    if is_sin == 1:
        return( np.sin(k*u) )
    else:
        return( np.cos(k*u) )

def get_theoretical_field( field, t ):
    amplitude = - epsilon * (m_e*c**2*k[field])/e * np.sin(wp*t)
    sin_flag = sin[field]
    x_contribution = get_contribution( sin_flag[0], kx )
    y_contribution = get_contribution( sin_flag[1], ky )
    z_contribution = get_contribution( sin_flag[2], kz )

    E = amplitude * x_contribution[:, np.newaxis, np.newaxis] \
                  * y_contribution[np.newaxis, :, np.newaxis] \
                  * z_contribution[np.newaxis, np.newaxis, :]

    return( E )


# Plot a single image of the fields at iteration 40
iteration = 40
ds = yt.load('./plt%05d/' %iteration)
t = ds.current_time.to_ndarray().mean()
data = ds.covering_grid(level=0, left_edge=ds.domain_left_edge,
                                dims=ds.domain_dimensions)
field ='Ez'
E_sim = data[field].to_ndarray()
plt.figure()
plt.imshow( E_sim[:,:,32] )
plt.colorbar()
plt.savefig('Plasma_wave_pattern.png')

# Get the time evolution of the plasma wave at a given point
Eval_sim = []
Eval_th = []
for iteration in tqdm(range(0,40,1)):
    ds = yt.load('./plt%05d/' %iteration)
    t = ds.current_time.to_ndarray().mean()
    data = ds.covering_grid(level=0, left_edge=ds.domain_left_edge,
                                        dims=ds.domain_dimensions)
    field ='Ez'
    E_sim = data[field].to_ndarray()
    E_th = get_theoretical_field(field, t)
    Eval_sim.append( E_sim[5,5,5] )
    Eval_th.append( E_th[5,5,5])

np.save('E_sim.npy', np.array(E_sim))
np.save('E_th.npy', np.array(E_th))

plt.clf()
plt.plot( Eval_sim )
plt.plot( Eval_th, 'r--')
plt.show()
