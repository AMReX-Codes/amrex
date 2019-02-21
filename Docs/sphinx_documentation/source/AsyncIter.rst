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
When the computations on a FAB finish, the ++ operator is called.
We overload this operator to traverse to the next runnable FAB.

Note: RGIter also supports data tiling.
Specifically, we overload the ++ operator so that it will traverse data tiles in a FAB before it goes to next FAB if the tiling flag in the FAB is enabled.
Instead of applying the computations in the loop body on the entire FAB, it executes them on a single tile at a time.


Generated Task Graph Code
=========================

The real input to the runtime system is an AMR program containing task dependency graphs (or task graph for short).
Thus, the code written with the above asynchronous iterators will be transformed into a task graph form.
The definition of a task dependency graph is as follows.
Each task of a graph performs some computations on an FArrayBox (FAB).
Tasks are connected with each other via edges, denoting the dependency on data.
A task can be executed when all data dependencies have been satisfied.
The code snippet below queries runnable tasks of a task dependency graph named regionGraph.
Note that each task dependency graph is more or less a wrapper of a MultiFab.
In this example, a task of regionGraph computes the body code of the while loop to update the associated FAB.
Each task of this graph receives data arrived at the runtime system and injects the data into the associated FAB.
After updating FAB, it lets the runtime know about the change.
The runtime system uses AMR domain knowledge to establish data dependencies among tasks, and thus it can answer which tasks are runnable and how to update neighbor FABs when a current FAB changes.

.. highlight:: c++

::

    while(!regionGraph->isGraphEmpty())
    {
        f = regionGraph->getAnyFireableRegion();
	multifabCopyPull(..., f, ...); //inject arrived dependent data into the fab, if any
        syncWorkerThreads();
	...//compute on the fab f of multifab associated with coarseRegionGraph
        syncWorkerThreads();
        multifabCopyPush(..., f, ...); //tell the runtime that data of Fab f changed
        regionGraph->finalizeRegion(f)
    }

The process of learning the domain knowledge is as follows.
At the beginning of the program, the runtime extracts the metadata needed for establishing data dependencies among tasks of the same graph or between two different graphs.
Every time the AMR grid hierarchy changes (i.e. when a few or all AMR levels regrid), the runtime re-extracts the metadata to correct the task dependency graphs.
Once the metadata extraction completes, the runtime system invokes the computation on AMR levels (e.g., timeStep, initTimeStep, and postTimeStep).

Known Limitations
=================

To realize enough task parallelism, the runtime system constructs a task dependency graph for the whole coarse time step and executes it asynchronously to the completion of the step.
As a result, any request to regrid an AMR level must be foreseen before the execution of a coarse time step.
If there is a regridding request during the graph execution, the runtime system simply ignores it.
In the future we may relax this constraint in the programming model.
However, such a support would come at a significant performance cost due to the required checkpointing and rollback activities.

