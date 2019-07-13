#! /usr/bin/env python

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
energyE = np.sum(scc.epsilon_0/2*(Ex**2+Ey**2+Ez**2))
energyB = np.sum(1./scc.mu_0/2*(Bx**2+By**2+Bz**2))
energy_end = energyE + energyB

Reflectivity = energy_end/energy_start
Reflectivity_theory = 1.8015e-06

print("Reflectivity", Reflectivity)
print("Reflectivity_theory", Reflectivity_theory)

assert( Reflectivity < 105./100 * Reflectivity_theory )
    
