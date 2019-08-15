#! /usr/bin/env python

'''
Analysis script of a WarpX simulation of rigid injection in a boosted frame.

A Gaussian electron beam starts from -5 microns, propagates rigidly up to
20 microns after which it expands due to emittance only (the focal position is
20 microns). The beam width is measured after ~50 microns, and compared with
the theory (with a 5% error allowed).

The simulation runs in a boosted frame, and the analysis is done in the lab
frame, i.e., on the back-transformed diagnostics.
'''

import sys, os, yt, glob
import numpy as np
import scipy.constants as scc
import read_raw_data
yt.funcs.mylog.setLevel(0)

# filename = sys.argv[1]

def get_particle_field(snapshot, species, field):
    fn = snapshot + '/' + species
    files = glob.glob(os.path.join(fn, field + '_*'))
    files.sort()
    all_data = np.array([])
    for f in files:
        data = np.fromfile(f)
        all_data = np.concatenate((all_data, data))
    return all_data

# Read data from back-transformed diagnostics
snapshot = './lab_frame_data/snapshot00001'
header   = './lab_frame_data/Header'
allrd, info = read_raw_data.read_lab_snapshot(snapshot, header)
z = np.mean( get_particle_field(snapshot, 'beam', 'z') )
w = np.std ( get_particle_field(snapshot, 'beam', 'x') )

# initial parameters
z0 = 20.e-6
w0 = 1.e-6
theta0 = np.arcsin(0.1)

# Theoretical beam width after propagation if rigid ON
wth = np.sqrt( w0**2 + (z-z0)**2*theta0**2 )
error = np.abs((w-wth)/wth)

# Print error and assert small error
print("Beam position: " + str(z))
print("Beam width   : " + str(w))
print("error: " + str(error))
assert( error < 0.03 )
