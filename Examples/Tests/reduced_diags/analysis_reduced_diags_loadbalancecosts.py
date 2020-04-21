#! /usr/bin/env python

# Copyright 2019-2020 Michael Rowan
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

# This script tests the reduced diagnostics `LoadBalanceCosts`.
# The setup is a uniform plasma with electrons.
# A measure of cost (based on timers or heuristic measure from cell and particles)
# diagnostic is output in the reduced diagnostic.  The efficiency (mean of cost
# per rank, normalized to the maximum cost over all ranks) extracted from the
# reduced diagnostic is compared before and after the load balance step; the test
# ensures that efficiency, measured via the reduced diagnostic, improves after
# the load balance step.

# Possible running time: ~ 1 s

import numpy as np
import sys

# Command line argument
fn = sys.argv[1]

# Load costs data
data = np.genfromtxt("./diags/reducedfiles/LBC.txt")
data = data[:,2:]

# Compute the number of datafields saved per box
n_data_fields = 0
with open("./diags/reducedfiles/LBC.txt") as f:
    h = f.readlines()[0]
    unique_headers=[''.join([l for l in w if not l.isdigit()]) for w in h.split()][2::]
    n_data_fields = len(set(unique_headers))
    f.close()

# From data header, data layout is:
#     [step, time,
#      cost_box_0, proc_box_0, lev_box_0, i_low_box_0, j_low_box_0, k_low_box_0(, gpu_ID_box_0 if GPU run), hostname_box_0,
#      cost_box_1, proc_box_1, lev_box_1, i_low_box_1, j_low_box_1, k_low_box_1(, gpu_ID_box_1 if GPU run), hostname_box_1,
#      ...
#      cost_box_n, proc_box_n, lev_box_n, i_low_box_n, j_low_box_n, k_low_box_n(, gpu_ID_box_n if GPU run), hostname_box_n]

# Function to get efficiency at an iteration i
def get_efficiency(i):
    # First get the unique ranks
    costs, ranks = data[i,0::n_data_fields], data[i,1::n_data_fields].astype(int)
    rank_to_cost_map = {r:0. for r in set(ranks)}

    # Compute efficiency before/after load balance and check it is improved
    for c, r in zip(costs, ranks):
        rank_to_cost_map[r] += c

    # Normalize the costs
    efficiencies = np.array(list(rank_to_cost_map.values()))
    efficiencies /= efficiencies.max()

    return efficiencies.mean()

# The iteration i=2 is load balanced; examine before/after load balance
efficiency_before, efficiency_after = get_efficiency(1), get_efficiency(2)
print('load balance efficiency (before load balance): ', efficiency_before)
print('load balance efficiency (after load balance): ', efficiency_after)

# The load balanced case is expcted to be more efficient then non-load balanced case
assert(efficiency_before < efficiency_after)
