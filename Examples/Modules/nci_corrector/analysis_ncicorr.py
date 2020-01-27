#! /usr/bin/env python

# Copyright 2019 Jean-Luc Vay, Maxence Thevenet, Remi Lehe
# Weiqun Zhang
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL


import sys
import yt
import re
import numpy as np
import scipy.constants as scc
yt.funcs.mylog.setLevel(0)

fn = sys.argv[1]
use_MR = re.search( 'nci_correctorMR', fn ) != None

if use_MR:
    energy_corrector_off = 5.e32
    energy_threshold = 1.e28
else:
    energy_corrector_off = 1.5e26
    energy_threshold = 1.e24

# Check EB energy after 1000 timesteps
filename = sys.argv[1]

ds = yt.load( filename )
ad0 = ds.covering_grid(level=0,left_edge=ds.domain_left_edge, dims=ds.domain_dimensions)
ex = ad0['boxlib', 'Ex'].v
ez = ad0['boxlib', 'Ez'].v
by = ad0['boxlib', 'By'].v
energy = np.sum(ex**2 + ez**2 + scc.c**2*by**2)

print("use_MR: %s" %use_MR)
print("energy if corrector off (from benchmark): %s" %energy_corrector_off)
print("energy threshold (from benchmark): %s" %energy_threshold)
print("energy from this run: %s" %energy)

assert( energy < energy_threshold )
