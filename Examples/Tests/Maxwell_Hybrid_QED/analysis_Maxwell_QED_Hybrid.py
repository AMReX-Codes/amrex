#!/usr/bin/env python3
chmod a+x analysis_Maxwell_QED_hyrbid.py
import sys
import yt ; ty.funcs.mylog.setLevel(0)
import numpy as np

QEDfile = './plt'
time = '00300'
dsQED = yt.load(QEDfile + time)
QED_all_data_level_0 = dsQED.covering_grid(level=0,left_edge=(dsQED.domain_left_edge),
                                           dims=dsQED.domain_dimensions)
EyQED = QED_all_data_level_0['boxlib', 'Ey'].v.squeeze()
resolution = dsQED.domain_dimensions[1]
dx = ((-0.5*resolution)+np.where(EyQED[0,:] == np.max(EyQED[0,:]))[0])*1024/resolution
v = dx/(float(dsQED.current_time))
vtheory = 299792458/np.sqrt((1+(12*10**(-13))/(8.8541878128*10**(-12)))/(1+(4*10**(-13))/(8.8541878128*10**(-12))))

print('Simulation velocity: ', v)
print('Theory veclority: ', vtheory)

assert(100*np.abs(v-vtheory)/vtheory < 1.25)

