#! /usr/bin/env python

import sys
import yt ; yt.funcs.mylog.setLevel(0)
import numpy as np
from scipy import signal

# Build Jx without filter (from other simulation)
my_F_nofilter = np.zeros([16,16])
my_F_nofilter[8,8] = -1.601068065642412e-11
my_F_nofilter[8,7] = -1.601068065642412e-11

# Build 2D filter
filter0 = np.array([.25,.5,.25])
my_order = [1,5]
my_filterx = filter0
my_filtery = filter0
while my_order[0]>1:
    my_filterx = np.convolve(my_filterx,filter0)
    my_order[0] -= 1
while my_order[1]>1:
    my_filtery = np.convolve(my_filtery,filter0)
    my_order[1] -= 1
my_filter = my_filterx[:,None]*my_filtery

# Apply filter. my_F_filetered is the theoretical value for filtered field
my_F_filtered = signal.convolve2d(my_F_nofilter, my_filter, boundary='symm', mode='same')

# Get simulation result for F_filtered
filename = sys.argv[1]
ds = yt.load( filename )
sl = yt.SlicePlot(ds, 2, 'jx', aspect=1)
all_data_level_0 = ds.covering_grid(level=0,left_edge=ds.domain_left_edge, dims=ds.domain_dimensions)
F_filtered = all_data_level_0['boxlib', 'jx'].v.squeeze()

# Compare theory and PIC for filtered value
error = np.sum( np.abs(F_filtered - my_F_filtered) ) / np.sum( np.abs(my_F_filtered) )
assert( error < 1.e-14 )
