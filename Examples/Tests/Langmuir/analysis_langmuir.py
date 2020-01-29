#! /usr/bin/env python

# Copyright 2019 Jean-Luc Vay, Maxence Thevenet, Remi Lehe
#
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL


import sys
import re
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.constants import c, e, m_e, epsilon_0
import numpy as np
import yt
yt.funcs.mylog.setLevel(50)

# this will be the name of the plot file
fn = sys.argv[1]
# Parse the direction in which the Langumir wave is launched
direction = re.search( 'Langmuir_([xyz])', fn ).group(1)

# Parameters of the plasma
u = 0.01
n0 = 1.e25
wp = (n0*e**2/(m_e*epsilon_0))**.5

# Load the dataset
ds = yt.load(fn)
t = ds.current_time.to_ndarray().mean() # in order to extract a single scalar
data = ds.covering_grid( 0, ds.domain_left_edge, ds.domain_dimensions )

# Check the J field along the direction of the wave, which oscillates at wp
j_predicted = -n0*e*c*u*np.cos( wp*t*39.5/40 ) # 40 timesteps / j at half-timestep
# Because of the shape factor, there are 2 cells that are incorrect
# at the edges of the plasma
if direction == 'x':
    j = data[ 'jx' ].to_ndarray()
    assert np.allclose( j[2:30,:,:], j_predicted, rtol=0.2 )
    assert np.allclose( j[34:-2,:,:], 0, atol=1.e-2 )
elif direction == 'y':
    j = data[ 'jy' ].to_ndarray()
    assert np.allclose( j[:,2:30,:], j_predicted, rtol=0.2 )
    assert np.allclose( j[:,34:-2,:], 0, atol=1.e-2 )
elif direction == 'z':
    j = data[ 'jz' ].to_ndarray()
    assert np.allclose( j[:,:,2:30], j_predicted, rtol=0.2 )
    assert np.allclose( j[:,:,34:-2], 0, atol=1.e-2 )

# Check the E field along the direction of the wave, which oscillates at wp
E_predicted = m_e * wp * u * c / e * np.sin(wp*t)
# Because of the shape factor, there are 2 cells that are incorrect
# at the edges of the plasma
if direction == 'x':
    E = data[ 'Ex' ].to_ndarray()
    # Print errors, and assert small error
    print( "relative error: np.max( np.abs( ( E[2:30,:,:] - E_predicted ) / E_predicted ) ) = %s" \
               %np.max( np.abs( ( E[2:30,:,:] - E_predicted ) / E_predicted ) ) )
    assert np.allclose( E[2:30,:,:], E_predicted, rtol=0.1 )
    print( "absolute error: np.max( np.abs( E[34:-2,:,:] ) ) = %s" %np.max( np.abs( E[34:-2,:,:] ) ) )
    assert np.allclose( E[34:-2,:,:], 0, atol=5.e-5 )
elif direction == 'y':
    E = data[ 'Ey' ].to_ndarray()
    # Print errors, and assert small error
    print( "relative error: np.max( np.abs( ( E[:,2:30,:] - E_predicted ) / E_predicted ) ) = %s" \
               %np.max( np.abs( ( E[:,2:30,:] - E_predicted ) / E_predicted ) ) )
    assert np.allclose( E[:,2:30,:], E_predicted, rtol=0.1 )
    print( "absolute error: np.max( np.abs( E[:,34:-2,:] ) ) = %s" %np.max( np.abs( E[:,34:-2,:] ) ) )
    assert np.allclose( E[:,34:-2,:], 0, atol=2.e-5 )
elif direction == 'z':
    E = data[ 'Ez' ].to_ndarray()
    # Print errors, and assert small error
    print( "relative error: np.max( np.abs( ( E[:,:,2:30] - E_predicted ) / E_predicted ) ) = %s" \
               %np.max( np.abs( ( E[:,:,2:30] - E_predicted ) / E_predicted ) ) )
    assert np.allclose( E[:,:,2:30], E_predicted, rtol=0.1 )
    print( "absolute error: np.max( np.abs( E[:,:,34:-2] ) ) = %s" %np.max( np.abs( E[:,:,34:-2] ) ) )
    assert np.allclose( E[:,:,34:-2], 0, atol=2.e-5 )

# Save an image to be displayed on the website
t_plot = np.linspace(0.0, t, 200)
plt.subplot(211)
plt.plot( t_plot, -n0 * e * c * u * np.cos( wp*t_plot ) )
plt.plot( 39.5/40*t, j_predicted, 'o' )
plt.ylabel( 'j'+direction )
plt.xlabel( 'Time' )
plt.subplot(212)
plt.plot( t_plot, m_e * wp * u * c / e * np.sin( wp*t_plot ) )
plt.plot( t, E_predicted, 'o' )
plt.ylabel( 'E'+direction )
plt.xlabel( 'Time' )
plt.savefig("langmuir_%s_analysis.png" %direction)
