.. sec:InputsTimeStepping:

Time Stepping
=============

The following inputs must be preceded by "amr."   Note that if both are specified, both criteria
are used and the simulation still stop when the first criterion is hit.  In the case of unsteady flow,
the simulation will stop when either the number of steps reaches max_step or time reaches stop_time.
In the case of unsteady flow, the simulation will stop when either the tolerance (difference between
subsequent steps) is reached or the number of iterations reaches the maximum number specified.

+------------------+-----------------------------------------------------------------------+-------------+-----------+
|                  | Description                                                           |   Type      | Default   |
+==================+=======================================================================+=============+===========+
| max_step         | Maximum number of time steps to take                                  |    Int      |  -1       |
+------------------+-----------------------------------------------------------------------+-------------+-----------+
| stop_time        | Maximum time to reach                                                 |    Real     | -1.0      |
+------------------+-----------------------------------------------------------------------+-------------+-----------+
