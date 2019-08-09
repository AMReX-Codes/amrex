#! /usr/bin/env python

import sys
import yt
import numpy as np
import scipy.constants as scc
yt.funcs.mylog.setLevel(0)

filename = sys.argv[1]

ds = yt.load( filename )
ad = ds.all_data()
ilev = ad['ions', 'particle_ionization_level'].v

N5_fraction = ilev[ilev == 5].size/ilev.size

print("Number of ions: " + str(ilev.size))
print("Number of N5+ : " + str(ilev[ilev == 5].size))
print("N5_fraction: " + str(N5_fraction))

assert ((N5_fraction > 0.30) and (N5_fraction < 0.34))
