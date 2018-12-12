.. _Chap:AsyncIter:

Amr Task
===============

``amrex/Src/AmrTask`` implements a runtime system that can execute the workload within each coarse time step in an asynchronous fashion, thereby overlapping communication with computation.
The runtime is defined in class RTS, which is a black box to the application developer.
Currently there are 3 implementations for RTS, located in ``amrex/Src/AmrTask/rts_impls``: MPI+Omp, MPI+Pthreads, and UPCXX+Pthreads.
The programmer can choose one of these backends to compile with their application without changing the application source code.
The following code snippet shows the main function of a typical AMReX application running asynchronously under the control of the runtime system.
The programmer first needs to use the namespace ``perilla``, which covers all the C++ classes for the runtime system.
To execute an AMR program (i.e. object of the Amr class), the programmer can simply create an object of RTS and pass the program object into the ``Iterate`` method. 
The runtime system will iteratively execute coarse time steps until the program completes.

.. highlight:: c++

::

    using namespace perilla;
    int main (int argc, char* argv[])
    {
        amrex::Initialize(argc,argv);
        ... //set up program input, e.g. start_time, stop_time, max_step
        Amr amr;
        amr.init(start_time,stop_time);
        RTS rts;
        rts.Iterate(&amr, max_step, stop_time);
        amrex::Finalize();
        return 0;
    }

In order to overlap communication with computation, the runtime system requires a sufficient amount of task parallelism.
Fortunately, there are plenty of opportunities to identify tasks in an AMReX program, including data parallel tasks (same workload on different data partitions) and control parallel tasks (different types of workload). 
To quickly create task parallelism that the runtime can exploit, the application can use asynchronous iterators that we develop in AMReX for interfacing between this AMR library and the runtime system. 
The API is very simple, and the asynchronous code is very similar to the original code using the synchronous multifab iterator (MFIter) described earlier in chapter Basics.

Before presenting the API of asynchronous iterators, we show how RTS executes generic task-based workload in an AMR coarse time step. 
The input to the runtime is an AMR program containing task graphs.
Each task of a graph performs some computations on an FArrayBox (FAB).
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
Note that the runtime exposes multiple threads per process, so the programmer needs to place sufficient memory protection for shared data within the process, e.g. when updating the state data.
In the code example above, syncWorkerThreads() synchronizes all threads exposed to the program.
This excludes internal threads that the runtime uses for handling communication.
This multithreaded interface adds some programming cost, but is necessary for mitigating the task scheduling overhead.
To avoid these programming details, the programmer can use of built-in iterators, such as fillpatch iterator and task graph iterator that we next discuss.

.. toctree::
   :maxdepth: 1


   AsyncIter
