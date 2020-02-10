#! /usr/bin/env python

# Copyright 2019-2020 Yinjian Zhao
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

# This script tests the reduced diagnostics.
# The setup is a uniform plasma with electrons and protons.
# Particle energy and field energy will be outputed
# using the reduced diagnostics.
# And they will be compared with the data in the plotfiles.

# Tolerance: 1.0e-8 for particle energy, 1.0e-3 for field energy.
# The difference of the field energy is relatively large,
# because fields data in plotfiles are cell-centered,
# but fields data in reduced diagnostics are staggered.

# Possible running time: ~ 2 s

import sys
import yt
import numpy as np
import scipy.constants as scc

fn = sys.argv[1]

ds = yt.load(fn)
ad = ds.all_data()

# PART1: get results from plotfiles

EPyt = 0.0
# electron
px = ad['electrons','particle_momentum_x'].to_ndarray()
py = ad['electrons','particle_momentum_y'].to_ndarray()
pz = ad['electrons','particle_momentum_z'].to_ndarray()
w  = ad['electrons','particle_weight'].to_ndarray()
EPyt = EPyt + np.sum( (np.sqrt((px**2+py**2+pz**2)*scc.c**2+scc.m_e**2*scc.c**4)-scc.m_e*scc.c**2)*w )
# proton
px = ad['protons','particle_momentum_x'].to_ndarray()
py = ad['protons','particle_momentum_y'].to_ndarray()
pz = ad['protons','particle_momentum_z'].to_ndarray()
w  = ad['protons','particle_weight'].to_ndarray()
EPyt = EPyt + np.sum( (np.sqrt((px**2+py**2+pz**2)*scc.c**2+scc.m_p**2*scc.c**4)-scc.m_p*scc.c**2)*w )

ad = ds.covering_grid(level=0, left_edge=ds.domain_left_edge, dims=ds.domain_dimensions)
Ex = ad['Ex'].to_ndarray()
Ey = ad['Ey'].to_ndarray()
Ez = ad['Ez'].to_ndarray()
Bx = ad['Bx'].to_ndarray()
By = ad['By'].to_ndarray()
Bz = ad['Bz'].to_ndarray()
Es = np.sum(Ex**2)+np.sum(Ey**2)+np.sum(Ez**2)
Bs = np.sum(Bx**2)+np.sum(By**2)+np.sum(Bz**2)
N  = np.array( ds.domain_width / ds.domain_dimensions )
dV = N[0]*N[1]*N[2]
EFyt = 0.5*Es*scc.epsilon_0*dV + 0.5*Bs/scc.mu_0*dV

# PART2: get results from reduced diagnostics

EFdata = np.genfromtxt("./diags/reducedfiles/EF.txt")
EPdata = np.genfromtxt("./diags/reducedfiles/EP.txt")
EF = EFdata[1][2]
EP = EPdata[1][2]

# PART3: print and assert

print('difference of field energy:', abs(EFyt-EF))
print('tolerance of field energy:', 1.0e-3)
print('difference of particle energy:', abs(EPyt-EP))
print('tolerance of particle energy:', 1.0e-8)

assert(abs(EFyt-EF) < 1.0e-3)
assert(abs(EPyt-EP) < 1.0e-8)
