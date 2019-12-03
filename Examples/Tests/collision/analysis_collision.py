#! /usr/bin/env python

# This script tests the collision module
# using electron-ion temperature relaxation.

# Possible errors:
# tolerance: 0.001
# Possible running time: ~ 30.0 s

import sys
import yt
import re
import math
import statistics

tolerance = 0.001

nt = 151
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

lastfilename = sys.argv[1]
temp = re.compile("([a-zA-Z_]+)([0-9]+)")
res = temp.match(lastfilename).groups()

j = 0
error = 0.0
for i in range (150,0,-1):
    # obtain file name
    if i > 140:
        filename = res[0]+"0000"+str(150-i)
    elif i > 50:
        filename = res[0]+"000" +str(150-i)
    else:
        filename = res[0]+"00"  +str(150-i)
    # load file
    ds  = yt.load( filename )
    ad  = ds.all_data()
    px  = ad['particle_momentum_x'].to_ndarray()
    # compute error
    vxe = statistics.mean(px[ 0:ne])/me/c
    vxi = statistics.mean(px[ne:np])/mi/c
    vxd = vxe - vxi
    fit = a*math.exp(b*j)
    error = error + abs(fit-vxd)
    j = j + 1

error = error / nt

print('error = ', error)
print('tolerance = ', tolerance)
assert(error < tolerance)
