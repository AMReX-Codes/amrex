# PrefixSum_MultiFab

This example implements a parallel prefix sum over MPI for 1-D MultiFab data.

The data to be summed is a field phi set to a uniform value of 1.0,
and the prefix sums at each point from the low end of the domain to the
high end are written out in a plotfile.

If `reverse=1` in the inputs, the prefix sum is done in order from
high to low ends of the domain instead.
