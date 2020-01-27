#! /usr/bin/env python

# Copyright 2019-2020 Luca Fedeli, Maxence Thevenet
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL


'''
Analysis script of a WarpX simulation of rigid injection.

A Gaussian electron beam starts from -5 microns, propagates rigidly up to
20 microns after which it expands due to emittance only (the focal position is
20 microns). The beam width is measured after ~50 microns, and compared with
the theory (with a 5% error allowed).

As a help to the user, the script also compares beam width to the theory in
case rigid injection is OFF (i.e., the beam starts expanding from -5 microns),
in which case a warning is raised.
'''

import sys
import yt
import numpy as np
yt.funcs.mylog.setLevel(0)

filename = sys.argv[1]

# WarpX headers include more data when rigid injection is used,
# which gives an error with the last yt release.
# To avoid this issue, the three last lines of WarpXHeader are removed if
# needed.
def remove_rigid_lines(plotfile, nlines_if_rigid):
    header_name = plotfile + '/WarpXHeader'
    f = open(header_name, 'r')
    file_lines = f.readlines()
    nlines = len(file_lines)
    f.close()
    if nlines == nlines_if_rigid:
        f = open(header_name, 'w')
        f.writelines(file_lines[:-3])
        f.close()

# Remove rigid injection header lines
remove_rigid_lines(filename, 18)
# Read beam parameters
ds = yt.load( filename )
ad = ds.all_data()
# Beam longitudinal position
z = np.mean(ad['beam', 'particle_position_y'].v)
# Beam width
w = np.std(ad['beam', 'particle_position_x'].v)

# initial parameters
z0 = 20.e-6
z0_no_rigid = -5.e-6
w0 = 1.e-6
theta0 = np.arcsin(0.1)

# Theoretical beam width after propagation if rigid OFF
# Inform the user if rigid injection simply off (just to be kind)
wth_no_rigid = np.sqrt( w0**2 + (z-z0_no_rigid)**2*theta0**2 )
error_no_rigid = np.abs((w-wth_no_rigid)/wth_no_rigid)
if ( error_no_rigid < 0.05):
    print("error no rigid: " + str(error_no_rigid))
    print("Looks like the beam defocuses as if rigid injection were OFF")

# Theoretical beam width after propagation if rigid ON
wth = np.sqrt( w0**2 + (z-z0)**2*theta0**2 )
error = np.abs((w-wth)/wth)
# Print error and assert small error
print("Beam position: " + str(z))
print("Beam width   : " + str(w))
print("error: " + str(error))
assert( error < 0.05 )
