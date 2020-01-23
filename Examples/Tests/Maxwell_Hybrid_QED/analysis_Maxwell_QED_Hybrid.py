#!/usr/bin/env python3
import sys
import yt ; yt.funcs.mylog.setLevel(0)
import numpy as np

dsQED = yt.load(sys.argv[1])
QED_all_data_level_0 = dsQED.covering_grid(level=0,left_edge=(dsQED.domain_left_edge),
                                           dims=dsQED.domain_dimensions)
EyQED = QED_all_data_level_0['boxlib', 'Ey'].v.squeeze()
resolution = dsQED.domain_dimensions[1]

scale = 1.e-6 # NOTE THIS IS DETERMINED FROM THE SIM SCALE SET IN INPUTS_2D. CHANGE ONE AND YOU MUST CHANGE THE OTHER

dx = scale*((-0.5*resolution)+np.where(EyQED[0,:] == np.max(EyQED[0,:]))[0])*1024/resolution
v = dx/(float(dsQED.current_time))
eps_0 = 8.8541878128*10**(-12)
xi = 10**(-23)

Es = 10**(5) #THIS IS THE STATIC FIELD FOUND IN INPUTS_2D. IF YOU CHANGE ONE AND YOU MUST CHANGE THE OTHER

c = 299792458.
vtheory = c/sqrt((1+12*xi*Es**2/eps_0)/(1+4*xi*Es**2/eps_0))

print('Simulation velocity: ', v)
print('Theory veclority: ', vtheory)

assert(100*np.abs(v-vtheory)/vtheory < 1.25)

