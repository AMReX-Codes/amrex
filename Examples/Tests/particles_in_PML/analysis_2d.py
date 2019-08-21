"""
This script tests the absorption of particles in the PML.

The input file inputs_2d is used: it features a positive and a
negative particle, going in opposite direction and eventually 
leaving the box. This script tests that the field in the box
is close to 0 once the particles have left. With regular 
PML, this test fails, since the particles leave a spurious
charge, with associated fields, behind them.
"""
import sys
import yt
import numpy as np
yt.funcs.mylog.setLevel(0)

# Open plotfile specified in command line
filename = sys.argv[1]
ds = yt.load( filename )

# Check that the field is low enough
ad0 = ds.covering_grid(level=0, left_edge=ds.domain_left_edge, dims=ds.domain_dimensions)
Ex_array = ad0['Ex'].to_ndarray()
Ey_array = ad0['Ey'].to_ndarray()
Ez_array = ad0['Ez'].to_ndarray()
max_Efield = max(Ex_array.max(), Ey_array.max(), Ez_array.max())
assert max_Efield < 0.0003