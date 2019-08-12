#! /usr/bin/env python

"""
This script tests the result of the ionization module in WarpX.

Input files inputs.rt and inputs.bf.rt are used to reproduce the test from
Chen, JCP, 2013, figure 2 (in the lab frame and in a boosted frame, 
respectively): a plane-wave laser pulse propagates through a
uniform N2+ neutral plasma and further ionizes the Nitrogen atoms. This test
checks that, after the laser went through the plasma, ~32% of Nitrogen
ions are N5+, in agreement with theory from Chen's article.
""" 

import sys
import yt
import numpy as np
import scipy.constants as scc
yt.funcs.mylog.setLevel(0)

# Open plotfile specified in command line, and get ion's ionization level.
filename = sys.argv[1]
ds = yt.load( filename )
ad = ds.all_data()
ilev = ad['ions', 'particle_ionization_level'].v

# Fraction of Nitrogen ions that are N5+.
N5_fraction = ilev[ilev == 5].size/float(ilev.size)

print("Number of ions: " + str(ilev.size))
print("Number of N5+ : " + str(ilev[ilev == 5].size))
print("N5_fraction   : " + str(N5_fraction))

assert ((N5_fraction > 0.30) and (N5_fraction < 0.34))
