#! /usr/bin/env python

# Copyright 2019 Andrew Myers, David Grote, Jean-Luc Vay
# Maxence Thevenet, Remi Lehe
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL


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

it = int(fn[-5:])

# Machine precision of the simulation
if ds.index._dtype == "float32":
    dfloat = 0.7e-4   # single: ok, just somewhat larger than 1.e-6 (6 digits)
else:
    dfloat = 2.0e-14  # double: ok, just a little larger than 1.e-15 (15 digits)

# Check the Jx field, which oscillates at wp
j_predicted = -n0*e*c*ux*np.cos( wp*t*(it-0.5)/it ) # j at half-timestep
jx = data['jx'].to_ndarray()
# Print errors, and assert small error
print( "relative error: np.max( np.abs( ( jx[:32,:,0] - j_predicted ) / j_predicted ) ) = %s" \
           %np.max( np.abs( ( jx[:32,:,0] - j_predicted ) / j_predicted ) ) )
assert np.allclose( jx[:32,:,0], j_predicted, rtol=0.1 )
print( "absolute error: np.max( np.abs( jx[32:,:,0] ) ) = %s" %np.max( np.abs( jx[32:,:,0] ) ) )
assert np.allclose( jx[32:,:,0], 0, atol=1.e-2 )

# Check the Ex field, which oscillates at wp
E_predicted = m_e * wp * ux * c / e * np.sin(wp*t)
Ex = data['Ex'].to_ndarray()
# Print errors, and assert small error
print( "relative error: np.max( np.abs( ( Ex[:32,:,0] - E_predicted ) / E_predicted ) ) = %s" \
           %np.max( np.abs( ( Ex[:32,:,0] - E_predicted ) / E_predicted ) ) )
assert np.allclose( Ex[:32,:,0], E_predicted, rtol=0.1 )
# predicted values should be zero here, but floating-point noise from E_predicted areas will sneak in
print( "absolute error: np.max( np.abs( Ex[32:,:,0] ) ) = %s" %np.max( np.abs( Ex[32:,:,0] ) ) )
float_noise_zero = np.max( np.abs( E_predicted * dfloat ) )
print( "float_noise_zero = %s"  % float_noise_zero )
assert np.allclose( Ex[32:,:,0], 0., atol=float_noise_zero )

# Save an image to be displayed on the website
t_plot = np.linspace(0.0, t, 200)
plt.subplot(211)
plt.plot( t_plot, -n0 * e * c * ux * np.cos( wp*t_plot ) )
plt.plot( (it-0.5)/it*t, j_predicted, 'o' )
plt.ylabel( 'jx' )
plt.xlabel( 'Time' )
plt.subplot(212)
plt.plot( t_plot, m_e * wp * ux * c / e * np.sin( wp*t_plot ) )
plt.plot( t, E_predicted, 'o' )
plt.ylabel( 'Ex' )
plt.xlabel( 'Time' )
plt.savefig("langmuir2d_analysis.png")
