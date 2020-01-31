#! /usr/bin/env python

# Copyright 2019 Maxence Thevenet, Revathi Jambunathan
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL


'''
Analysis script of a WarpX simulation in a boosted frame.

The simulation runs in a boosted frame, and the analysis is done in the lab
frame, i.e., on the back-transformed diagnostics for the full 3D simulation and
an x-z slice at y=y_center. The field-data, Ez, along z, at (x_center,y_center,:) is compared
between the full back-transformed diagnostic and the reduced diagnostic (i.e., x-z slice) .
'''

import numpy as np
import read_raw_data

# Read data from back-transformed diagnostics of entire domain
snapshot = './lab_frame_data/snapshots/snapshot00002'
header   = './lab_frame_data/snapshots/Header'
allrd, info = read_raw_data.read_lab_snapshot(snapshot, header)
F = allrd['Ez']
print("F.shape ", F.shape)
F_1D = np.squeeze(F[F.shape[0]//2,F.shape[1]//2,:])


# Read data from reduced back-transformed diagnostics (i.e. slice)
snapshot_slice = './lab_frame_data/slices/slice00002'
header_slice   = './lab_frame_data/slices/Header'
allrd, info = read_raw_data.read_lab_snapshot(snapshot_slice, header_slice)
Fs = allrd['Ez']
print("Fs.shape", Fs.shape)
Fs_1D = np.squeeze(Fs[Fs.shape[0]//2,1,:])

error = np.max(np.abs(Fs_1D - F_1D)) / np.max(np.abs(F_1D))

# Print error and assert small error
print("relative error: " + str(error))
assert( error < 1E-15 )
