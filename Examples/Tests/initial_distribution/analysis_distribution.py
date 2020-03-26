#! /usr/bin/env python

# Copyright 2019-2020 Yinjian Zhao
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

# This script tests initial distributions.
# 1 denotes gaussian distribution.
# 2 denotes maxwell-boltzmann distribution.
# 3 denotes maxwell-juttner distribution.
# 4 denotes gaussian position distribution.
# The distribution is obtained through reduced diagnostic ParticleHistogram.

import numpy as np
import scipy.constants as scc
import scipy.special as scs
from read_raw_data import read_reduced_diags_histogram

# print tolerance
tolerance = 0.007
print('Tolerance:', tolerance)

#===============================
# gaussian and maxwell-boltzmann
#===============================

# load data
bin_value, h1x = read_reduced_diags_histogram("h1x.txt")[2:]
h1y = read_reduced_diags_histogram("h1y.txt")[3]
h1z = read_reduced_diags_histogram("h1z.txt")[3]
h2x = read_reduced_diags_histogram("h2x.txt")[3]
h2y = read_reduced_diags_histogram("h2y.txt")[3]
h2z = read_reduced_diags_histogram("h2z.txt")[3]

# parameters of theory
u_rms = 0.01
gamma = np.sqrt(1.0+u_rms*u_rms)
v_rms = u_rms / gamma * scc.c
n     = 1.0e21
V     = 8.0
db    = 0.0016

# compute the analytical solution
f = n*V*scc.c*db*np.exp(-0.5*(bin_value*scc.c/v_rms)**2)/(v_rms*np.sqrt(2.0*scc.pi))
f_peak = np.amax(f)

# compute error
# note that parameters are chosen such that gaussian and
# maxwell-boltzmann distributions are identical
f1_error = np.sum(np.abs(f-h1x)+np.abs(f-h1y)+np.abs(f-h1z))/bin_value.size / f_peak
f2_error = np.sum(np.abs(f-h2x)+np.abs(f-h2y)+np.abs(f-h2z))/bin_value.size / f_peak

print('Gaussian distribution difference:', f1_error)
print('Maxwell-Boltzmann distribution difference:', f2_error)

assert(f1_error < tolerance)
assert(f2_error < tolerance)

#================
# maxwell-juttner
#================

# load data
bin_value, bin_data = read_reduced_diags_histogram("h3.txt")[2:]

# parameters of theory
theta = 1.0
K2    = scs.kn(2,1.0/theta)
n     = 1.0e21
V     = 8.0
db    = 0.22

# compute the analytical solution

f = n*V*db * bin_value**2 * np.sqrt(1.0-1.0/bin_value**2) / \
    (theta*K2) * np.exp(-bin_value/theta)
f_peak = np.amax(f)

# compute error
f3_error = np.sum( np.abs(f-bin_data) ) / bin_value.size / f_peak

print('Maxwell-Juttner distribution difference:', f3_error)

assert(f3_error < tolerance)

#==============
# gaussian beam
#==============

# load data
bin_value, h4x = read_reduced_diags_histogram("h4x.txt")[2:]
h4y = read_reduced_diags_histogram("h4y.txt")[3]
h4z = read_reduced_diags_histogram("h4z.txt")[3]

# parameters of theory
x_rms = 0.25
q_tot = -1.0e-20
q_e   = -1.602176634e-19
npart = q_tot/q_e
db    = 0.04

# compute the analytical solution
f = npart*db * np.exp(-0.5*(bin_value/x_rms)**2)/(x_rms*np.sqrt(2.0*scc.pi))
f_peak = np.amax(f)

# compute error
f4_error = np.sum(np.abs(f-h4x)+np.abs(f-h4y)+np.abs(f-h4z))/bin_value.size / f_peak

print('Gaussian position distribution difference:', f4_error)

assert(f4_error < tolerance)
