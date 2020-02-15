#! /usr/bin/env python

# Copyright 2019 Andrew Myers, Jean-Luc Vay, Maxence Thevenet
# Remi Lehe, Weiqun Zhang
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL


import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

# this will be the name of the plot file
fn = sys.argv[1]

# you can save an image to be displayed on the website
t = np.arange(0.0, 2.0, 0.01)
s = 1 + np.sin(2*np.pi*t)
plt.plot(t, s)
plt.savefig("laser_analysis.png")

# return '0' for success, anything else for failure
sys.exit(0)
