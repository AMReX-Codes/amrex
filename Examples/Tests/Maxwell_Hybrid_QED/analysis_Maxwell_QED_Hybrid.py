#!/usr/bin/env python3

# Copyright 2019-2020 Axel Huebl, Glenn Richardson, Maxence Thevenet
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

import sys
import yt ; yt.funcs.mylog.setLevel(0)
import numpy as np
import scipy.constants as scc

# Static electric field and quantum parameters, from the input file.
Es = 1.0e5
xi = 1.0e-23

# Load dataset and get laser field
dsQED = yt.load(sys.argv[1])
QED_all_data_level_0 = dsQED.covering_grid(level=0,left_edge=(dsQED.domain_left_edge),
                                           dims=dsQED.domain_dimensions)
EyQED_2d = QED_all_data_level_0['boxlib', 'Ey'].v.squeeze()

# Extract 1D lineout of the laser field
EyQED = EyQED_2d[EyQED_2d.shape[0]//2,:]

# Longitudinal resolution
dz = dsQED.domain_width[1].v/dsQED.domain_dimensions[1]

# Initial position of the laser pulse max (from input file)
z_start = 0.
# Final position of the laser pulse max (from plotfile)
z_end = dsQED.domain_left_edge[1].v + np.argmax(EyQED) * dz
# Compute phase velocity and compare with theory
phase_velocity_pic = (z_end-z_start)/dsQED.current_time.v
phase_velocity_theory = scc.c/np.sqrt((1.+12.*xi*Es**2/scc.epsilon_0)/(1.+4.*xi*Es**2/scc.epsilon_0))
error_percent = 100.*np.abs(phase_velocity_pic-phase_velocity_theory)/phase_velocity_theory

# Print and assert correctness
print('Simulation velocity: ' + str(phase_velocity_pic))
print('Theory velocity    : ' + str(phase_velocity_theory))
print('error (%)          : ' + str(error_percent) )
print('Theoretical difference between with/without QED (%): ' + str(100*np.abs(phase_velocity_theory-scc.c)/scc.c))
assert( error_percent < 1.25 )
