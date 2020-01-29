#!/usr/bin/env python

# Copyright 2019-2020 Axel Huebl, Remi Lehe
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

"""
This script checks the space-charge initialization routine, by
verifying that the space-charge field of a Gaussian beam corresponds to
the expected theoretical field.
"""
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import yt
import numpy as np
import scipy.constants as scc
from scipy.special import gammainc
yt.funcs.mylog.setLevel(0)

# Parameters from the Simulation
Qtot = -1.e-20
r0 = 2.e-6

# Open data file
filename = sys.argv[1]
ds = yt.load( filename )
# Extract data
ad0 = ds.covering_grid(level=0, left_edge=ds.domain_left_edge, dims=ds.domain_dimensions)
Ex_array = ad0['Ex'].to_ndarray().squeeze()
if ds.dimensionality == 2:
    # Rename the z dimension as y, so as to make this script work for 2d and 3d
    Ey_array = ad0['Ez'].to_ndarray().squeeze()
elif ds.dimensionality == 3:
    Ey_array = ad0['Ey'].to_ndarray()
    Ez_array = ad0['Ez'].to_ndarray()

# Extract grid coordinates
Nx, Ny, Nz =  ds.domain_dimensions
xmin, ymin, zmin = ds.domain_left_edge.v
Lx, Ly, Lz = ds.domain_width.v
x = xmin + Lx/Nx*(0.5+np.arange(Nx))
y = ymin + Ly/Ny*(0.5+np.arange(Ny))
z = zmin + Lz/Nz*(0.5+np.arange(Nz))

# Compute theoretical field
if ds.dimensionality == 2:
    x_2d, y_2d = np.meshgrid(x, y, indexing='ij')
    r2 = x_2d**2 + y_2d**2
    factor = (Qtot/r0)/(2*np.pi*scc.epsilon_0*r2) * (1-np.exp(-r2/(2*r0**2)))
    Ex_th = x_2d * factor
    Ey_th = y_2d * factor
elif ds.dimensionality == 3:
    x_2d, y_2d, z_2d = np.meshgrid(x, y, z, indexing='ij')
    r2 = x_2d**2 + y_2d**2 + z_2d**2
    factor = Qtot/(4*np.pi*scc.epsilon_0*r2**1.5) * gammainc(3./2, r2/(2.*r0**2))
    Ex_th = factor*x_2d
    Ey_th = factor*y_2d
    Ez_th = factor*z_2d

# Plot theory and data
def make_2d(arr):
    if arr.ndim == 3:
        return arr[:,:,Nz//2]
    else:
        return arr
plt.figure(figsize=(10,10))
plt.subplot(221)
plt.title('Ex: Theory')
plt.imshow(make_2d(Ex_th))
plt.colorbar()
plt.subplot(222)
plt.title('Ex: Simulation')
plt.imshow(make_2d(Ex_array))
plt.colorbar()
plt.subplot(223)
plt.title('Ey: Theory')
plt.imshow(make_2d(Ey_th))
plt.colorbar()
plt.subplot(224)
plt.title('Ey: Simulation')
plt.imshow(make_2d(Ey_array))
plt.colorbar()
plt.savefig('Comparison.png')

# Automatically check the results
def check(E, E_th, label):
    print( 'Relative error in %s: %.3f'%(
            label, abs(E-E_th).max()/E_th.max()))
    assert np.allclose( E, E_th, atol=0.1*E_th.max() )

check( Ex_array, Ex_th, 'Ex' )
check( Ey_array, Ey_th, 'Ey' )
if ds.dimensionality == 3:
    check( Ez_array, Ez_th, 'Ez' )
