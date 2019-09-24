#! /usr/bin/env python

import sys
import yt
import numpy as np
import scipy.constants as scc
yt.funcs.mylog.setLevel(0)

# Check EB energy after 1000 timesteps
filename = sys.argv[1]
ds = yt.load( filename )
ad = ds.all_data()
ex = np.reshape(ad['boxlib', 'Ex'].v,(128,128))
ez = np.reshape(ad['boxlib', 'Ez'].v,(128,128))
by = np.reshape(ad['boxlib', 'By'].v,(128,128))
energy = np.sum(ex**2 + ez**2 + scc.c**2*by**2)*1.e-12

print("energy: %s" %energy)
print("Should be ~1.e0 if filter ON, ~1.e5 if filter OFF.")

assert( energy < 7. )
