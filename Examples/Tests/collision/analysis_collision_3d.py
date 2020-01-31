#! /usr/bin/env python

# Copyright 2019-2020 Yin-YinjiaZhao, Yinjian Zhao
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

# This script tests the collision module
# using electron-ion temperature relaxation in 3D.
# Initially, electrons and ions are both in equilibrium
# (gaussian) distributions, but have different temperatures.
# Relaxation occurs to bring the two temeratures to be
# a final same temperature through collisions.
# The code was tested to be valid, more detailed results
# were used to obtian an exponential fit with
# coefficients a and b.
# This automated test compares the results with the fit.

# Possible errors:
# tolerance: 0.001
# Possible running time: ~ 30.0 s

import sys
import yt
import re
import math
import numpy
from glob import glob

tolerance = 0.001

ng = 512
ne = ng * 200
ni = ng * 200
np = ne + ni

c  = 299792458.0
me = 9.10938356e-31
mi = me * 5.0

# exponential fit coefficients
a =  0.041817463099883
b = -0.083851393560288

last_fn = sys.argv[1]
temp = re.compile("([a-zA-Z_]+)([0-9]+)")
res = temp.match(last_fn).groups()
fn_list = glob(res[0] + "?????")

error = 0.0
nt = 0
for fn in fn_list:
    # load file
    ds  = yt.load( fn )
    ad  = ds.all_data()
    px  = ad['particle_momentum_x'].to_ndarray()
    # get time index j
    buf = temp.match(fn).groups()
    j = int(buf[1])
    # compute error
    vxe = numpy.mean(px[ 0:ne])/me/c
    vxi = numpy.mean(px[ne:np])/mi/c
    vxd = vxe - vxi
    fit = a*math.exp(b*j)
    error = error + abs(fit-vxd)
    nt = nt + 1

error = error / nt

print('error = ', error)
print('tolerance = ', tolerance)
assert(error < tolerance)
