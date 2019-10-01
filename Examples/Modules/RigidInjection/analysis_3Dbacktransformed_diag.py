#! /usr/bin/env python

'''
Analysis script of a WarpX simulation in a boosted frame.

The simulation runs in a boosted frame, and the analysis is done in the lab
frame, i.e., on the back-transformed diagnostics for the full 3D simulation and
an x-z slice at y=y_center. The field-data, Ez, along z, at (x_center,y_center,:) is compared
between the full back-transformed diagnostic and the reduced diagnostic (i.e., x-z slice) .
'''

import sys, os, yt, glob
import numpy as np
import scipy.constants as scc
import read_raw_data
yt.funcs.mylog.setLevel(0)

# Read data from back-transformed diagnostics of entire domain
snapshot = './lab_frame_data/snapshots/snapshot00000'
header   = './lab_frame_data/snapshots/Header'
allrd, info = read_raw_data.read_lab_snapshot(snapshot, header)
F = allrd['Ez']
F_1D = np.squeeze(F[F.shape[0]//2,F.shape[1]//2,:])


# Read data from reduced back-transformed diagnostics (i.e. slice)
snapshot_slice = './lab_frame_data/slices/slice00000'
header_slice   = './lab_frame_data/slices/Header'
allrd, info = read_raw_data.read_lab_snapshot(snapshot_slice, header_slice)
Fs = allrd['Ez']
Fs_1D = np.squeeze(Fs[Fs.shape[0]//2,Fs.shape[1]//2,:])

error = np.max(np.abs(Fs_1D - F_1D))

# Print error and assert small error
print("absolute error: " + str(error))
assert( error < 1E-15 )
