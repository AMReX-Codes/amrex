#! /usr/bin/env python
"""
This script tests the result of the Galilen method in WarpX.
It compares the energy of the electric field calculated using Galilean method with
'v_galiean = (0.,0., 0.99498743710662)' versus standard PSATD (v_galiean = (0.,0.,0.)):
    * if 'v_galilean == 0': simulation is unstable because of the arosen NCI;
    * if 'v_galilean != 0 : NCI is suppresed => simulation is stable.
"""
import sys
import yt ; yt.funcs.mylog.setLevel(0)
import numpy as np
import scipy.constants as scc

filename = sys.argv[1]


ds = yt.load( filename )
all_data_level_0 = ds.covering_grid(level=0,left_edge=ds.domain_left_edge, dims=ds.domain_dimensions)
Ex = all_data_level_0['boxlib', 'Ex'].v.squeeze()
Ey = all_data_level_0['boxlib', 'Ey'].v.squeeze()
Ez = all_data_level_0['boxlib', 'Ez'].v.squeeze()

#E field energy calculated with Galilean method (v_galilean = (0,0,0.99498743710662))
energyE_gal_psatd = np.sum(scc.epsilon_0/2*(Ex**2+Ey**2+Ez**2))

#E field energy precalculated with standard PSATD (v_galilean = (0,0,0))
energyE_psatd = 270975.396667626 #E field energy calculated with PSATD (v_galilean = (0,0,0))

assert( energyE_gal_psatd < 1e-10 * energyE_psatd )
