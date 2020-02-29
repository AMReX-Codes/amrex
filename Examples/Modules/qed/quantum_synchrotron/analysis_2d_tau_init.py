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

# This script checks if electrons and positrons initialized with
# Quantum Synchrotron process enabled
# do actually have an exponentially distributed optical depth

# Tolerance
tol = 1e-2

def check():
    filename = sys.argv[1]
    data_set = yt.load(filename)

    all_data = data_set.all_data()
    res_ele_tau = all_data["electrons", 'particle_optical_depth_QSR']
    res_pos_tau = all_data["positrons", 'particle_optical_depth_QSR']

    loc_ele, scale_ele = st.expon.fit(res_ele_tau)
    loc_pos, scale_pos = st.expon.fit(res_pos_tau)

    # loc should be very close to 0, scale should be very close to 1
    assert(np.abs(loc_ele - 0) < tol)
    assert(np.abs(loc_pos - 0) < tol)
    assert(np.abs(scale_ele - 1) < tol)
    assert(np.abs(scale_pos - 1) < tol)

def main():
    check()

if __name__ == "__main__":
    main()

