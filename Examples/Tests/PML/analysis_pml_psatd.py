#! /usr/bin/env python

# Copyright 2019 Jean-Luc Vay, Maxence Thevenet, Remi Lehe
#
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL


import sys
import yt ; yt.funcs.mylog.setLevel(0)
import numpy as np
import scipy.constants as scc

filename = sys.argv[1]

############################
### INITIAL LASER ENERGY ###
############################
energy_start = 9.1301289517e-08

##########################
### FINAL LASER ENERGY ###
##########################
ds = yt.load( filename )
all_data_level_0 = ds.covering_grid(level=0,left_edge=ds.domain_left_edge, dims=ds.domain_dimensions)
Bx = all_data_level_0['boxlib', 'Bx'].v.squeeze()
By = all_data_level_0['boxlib', 'By'].v.squeeze()
Bz = all_data_level_0['boxlib', 'Bz'].v.squeeze()
Ex = all_data_level_0['boxlib', 'Ex'].v.squeeze()
Ey = all_data_level_0['boxlib', 'Ey'].v.squeeze()
Ez = all_data_level_0['boxlib', 'Ez'].v.squeeze()
rho = all_data_level_0['boxlib','rho'].v.squeeze()
divE = all_data_level_0['boxlib','divE'].v.squeeze()
energyE = np.sum(scc.epsilon_0/2*(Ex**2+Ey**2+Ez**2))
energyB = np.sum(1./scc.mu_0/2*(Bx**2+By**2+Bz**2))
energy_end = energyE + energyB

Reflectivity = energy_end/energy_start
Reflectivity_theory = 1.3806831258153887e-06

print("Reflectivity: %s" %Reflectivity)
print("Reflectivity_theory: %s" %Reflectivity_theory)

assert( abs(Reflectivity-Reflectivity_theory) < 5./100 * Reflectivity_theory )

# Check relative L-infinity spatial norm of rho/epsilon_0 - div(E)
Linf_norm = np.amax( np.abs( rho/scc.epsilon_0 - divE ) ) / np.amax( np.abs( rho/scc.epsilon_0 ) )
assert( Linf_norm < 2.e-2 )
