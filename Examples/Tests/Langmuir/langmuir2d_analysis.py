#! /usr/bin/env python

import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.constants import c, e, m_e, epsilon_0
import numpy as np
import yt
yt.funcs.mylog.setLevel(50)

# this will be the name of the plot file
fn = sys.argv[1]

# Parameters of the plasma
ux = 0.01
n0 = 1.e25
wp = (n0*e**2/(m_e*epsilon_0))**.5

# Load the dataset
ds = yt.load(fn)
t = ds.current_time.to_ndarray().mean() # in order to extract a single scalar
data = ds.covering_grid( 0, ds.domain_left_edge, ds.domain_dimensions )

# Check the J fields
assert np.allclose( data['jz'].to_ndarray(), 0, atol=0.1 )
assert np.all( data['jy'].to_ndarray() == 0. )
# Check the Jx field, which oscillates at wp
j_predicted = -n0*e*c*ux*np.cos( wp*t*39.5/40 ) # 40 timesteps / j at half-timestep
jx = data['jx'].to_ndarray()
assert np.allclose( jx[:32,:,0], j_predicted, rtol=0.1 )
assert np.allclose( jx[32:,:,0], 0, atol=1.e-2 )

# Check the E fields
assert np.allclose( data['Ez'].to_ndarray(), 0, atol=1.e-4 )
assert np.all( data['Ey'].to_ndarray() == 0. )
# Check the Ex field, which oscillates at wp
E_predicted = m_e * wp * ux * c / e * np.sin(wp*t)
Ex = data['Ex'].to_ndarray()
assert np.allclose( Ex[:32,:,0], E_predicted, rtol=0.1 )
assert np.allclose( Ex[32:,:,0], 0, atol=1.e-4 )

# Check the B fields
assert np.all( data['Bx'].to_ndarray() == 0. )
assert np.allclose( data['By'].to_ndarray(), 0, rtol=2. )
assert np.all( data['Bz'].to_ndarray() == 0. )

# Save an image to be displayed on the website
t_plot = np.linspace(0.0, t, 200)
plt.subplot(211)
plt.plot( t_plot, -n0 * e * c * ux * np.cos( wp*t_plot ) )
plt.plot( 39.5/40*t, j_predicted, 'o' )
plt.ylabel( 'jx' )
plt.xlabel( 'Time' )
plt.subplot(212)
plt.plot( t_plot, m_e * wp * ux * c / e * np.sin( wp*t_plot ) )
plt.plot( t, E_predicted, 'o' )
plt.ylabel( 'Ex' )
plt.xlabel( 'Time' )
plt.savefig("langmuir2d_analysis.png")
