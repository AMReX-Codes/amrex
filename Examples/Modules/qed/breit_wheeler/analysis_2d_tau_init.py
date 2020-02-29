#! /usr/bin/env python

# Copyright 2019 Luca Fedeli, Maxence Thevenet
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

import yt
import numpy as np
import scipy.stats as st
import sys

# This script checks if photons initialized with Breit Wheeler process enabled
# do actually have an exponentially distributed optical depth

# Tolerance
tol = 1e-2

def check():
    filename = sys.argv[1]
    data_set = yt.load(filename)

    all_data = data_set.all_data()
    res_tau = all_data["photons", 'particle_optical_depth_BW']

    loc, scale = st.expon.fit(res_tau)

    # loc should be very close to 0, scale should be very close to 1
    assert(np.abs(loc - 0) < tol)
    assert(np.abs(scale - 1) < tol)

def main():
    check()

if __name__ == "__main__":
    main()

