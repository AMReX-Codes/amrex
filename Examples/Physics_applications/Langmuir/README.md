Simple example of Langmuir oscillations in a uniform plasma.

To see an example of how to use the level 0 dynamic load balancing, follow the following steps
in this directory.

1) Build an MPI-enabled executable.

2) Run with 4 processes using inputs.lb
   mv plt00010 to plt00010.lb

3) Run with 4 processes using inputs.nolb
   mv plt00010 to plt00010.nolb

4) amrvis3d plt00010.lb plt00010.nolb
   set the field to "part_per_proc"

You should see the effect of load balancing based of the number of particles.

In the nolb case, there are 4484 particles per process on two of the processes, and 644 on the other two.

In the   lb case, there are 2916, 2852, 2244, and 2244 particles per process.

