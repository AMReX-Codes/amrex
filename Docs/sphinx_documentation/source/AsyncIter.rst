FillPatch Iterator
==================

FillPatch is an important operation commonly used in AMReX applications.
This operation interpolates data in both space and time.
Communication between AMR levels may incur when FillPatch interpolates data from a coarse AMR level and stores the result on the next finer level.
This operation also results in communication within the same AMR level when the subcycling option is used, which requires data interpolation in time.

We develop an asynchronous version of the FillPatch operation, called Asynchronous FillPatch Iterator.
Each iterator takes care of the communication with the previous and next subcycles at the same AMR level (time) and between the current and the next finer AMR levels (space and time).
The iterator first automatically prepares temporary data needed for these communication activities and the data connections (aka data paths or data dependencies) among them.

Based on this setup, the programmer can design numerical solvers.
This work is fairly simple.
At a certain simulation time on an AMR level, the programmer can ask the runtime system which FABs have received sufficient data for advancing to the next time step.
Although the FillPatch operation can be handled independently by the communication handler of the runtime system, this operation requires some computations such as packing/unpacking and extrapolation.
The programmer has the freedom to dedicate a few threads from the pool of worker threads to parallelize those computations.
This design choice may help the runtime process FillPatch operations faster, but may slow down the main computation.
Thus, our advise to the programmer on using how many threads for the FillPatch is that it depends on the compute intensity of the actual workload.
If the simulation is memory-bandwidth or network-bandwidth bound, the programmer can get the benefit from sparing more threads for doing FillPatch.

RegionGraph Iterator
====================

We can simplify the programming work further with a new abstraction called RegionGraph Iterator a.k.a RGIter.
This abstraction is a for loop (see the following code snippet), which can hide details of the asynchronous FillPatch Iterator in the init part and the graph traversal in the ++ operator.
The only job required from the programmer is to specify the computations on the data, and they can easily place these computations in the loop body.

.. highlight:: c++

::

    for (RGIter rgi(afpi_vec, upper_afpi_vec, ...); rgi.isValid(); ++rgi){
        int f = rgi.currentRegion;
	...//computation on FAB f
    }

The execution of RGIter is as follows.
Initially, an object of RGIter (i.e. rgi) is instantiated, taking vectors of FillPatch Iterators on the current and upper AMR levels as arguments (each element of the vector corresponds to a specific time).
Based on these arguments, a task dependency graph spanning two AMR levels will be established. 
Next, isValid() asks the runtime system for FABs that have received all dependent data.
When there is such a FAB, the computations in the loop body can execute on the FAB's data.
When the computations on a FAB finishes, the ++ operator is called.
We overload this operator to traverse to the next runnable FAB.

Note: RGIter also supports data tiling.
Specifically, we overload the ++ operator so that it will traverse data tiles in a FAB before it goes to next FAB if the tiling flag in the FAB is enabled.
Instead of applying the computations in the loop body on the entire FAB, it executes them on a single tile at a time.

