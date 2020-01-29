#! /usr/bin/env python

# Copyright 2019 Yinjian Zhao
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL


# This script tests the particle pusher (HC)
# using a force-free field,
# in which position x should remain 0.
# An initial velocity Vy corresponding to
# Lorentz factor = 20 is used.
# Bz is fixed at 1 T.
# Ex = -Vy*Bz.

# Possible errors:
# Boris:     2321.3958529
# Vay:       0.00010467
# HC:        0.00011403
# tolerance: 0.001
# Possible running time: ~ 4.0 s

import sys
import yt

tolerance = 0.001

filename = sys.argv[1]
ds = yt.load( filename )
ad = ds.all_data()
x  = ad['particle_position_x'].to_ndarray()

print('error = ', abs(x))
print('tolerance = ', tolerance)
assert(abs(x) < tolerance)
