.. _Chap:AmrTask:

Amr Task
===============

``amrex/Src/AmrTask`` implements a runtime system that can execute the workload within each coarse time step in an asynchronous fashion, thereby overlapping communication with computation.
The runtime is defined in class RTS, which is a black box to the application developer.
Currently there are 3 implementations for RTS, located in ``amrex/Src/AmrTask/rts_impls``: MPI+Omp, MPI+Pthreads, and UPCXX+Pthreads.
The programmer can choose one of these backends to compile with their application without changing the application source code.
The following code snippet shows the main function of a typical AMR application running asynchronously under the control of the runtime system.
The programmer first needs to use the namespace :cpp:`perilla`, which covers all the C++ classes for the runtime system.
To execute an AMR program (i.e. object of the Amr class), the programmer can simply create an object of RTS and pass the program object into the :cpp:`Iterate` method. 
The runtime system will iteratively execute coarse time steps until the program completes.
In order to overlap communication with computation, the runtime requires enough task parallelism.
To quickly create task-level parallelism that the runtime can exploit, the application can use asynchronous grid iterators that we develop in AMReX. 
The API is very simple, and the asynchronous code is very similar to the original code using the synchronous multifab iterator (MFIter) previously mentioned in chapter Basics.
We begin with showing how RTS executes functions in an Amr coarse time step. 

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

The runtime system schedules AMR workload in a coarse time step as follows.
At the beginning of the program, the runtime extracts the metadata needed for establising data dependencies of task dependency graphs represented the application.
Everytime the AMR grid hierarchy changes when a few or all AMR levels regrid, the runtime re-extracts the metadata to correct the task depdenendy graphs.
Once the metadata extraction completes, the runtime system invokes the computation on AMR levels (e.g., timeStep, initTimeStep, and postTimeStep).
Note that the runtime exposes multiple threads per process, so the programmer needs to place sufficient memory protection for shared data within the process, e.g. when updating the state data.
This interface adds some programming cost, but is neccessary for mitigating the task scheduling overhead.
To avoid programming details, the programmer can make use of builtin iterators, such as fillpatch iterator task graph iterator.


Asynchronous fillpatch iterator
:cpp:`AsyncFillPatchIterator`


Region Graph iterator
:cpp:`RGIter`

.. toctree::
   :maxdepth: 1


   AmrTask
