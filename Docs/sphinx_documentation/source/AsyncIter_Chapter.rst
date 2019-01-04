.. _Chap:AsyncIter:

Asynchronous Iterators (AmrTask)
================================

Hiding communication overheads via overlapping communication with computation requires a sufficiently large amount of task parallelism.
This problem is even more challenging due to various types of tasks in an AMReX program, including data parallel tasks (same workload on different data partitions) and control parallel tasks (different types of workload). 
This chapter introduces the API of AMReX's asynchronous iterators that can facilitate the job of identifying tasks in the applications.  
We have developed two iterators called FillPatch and RegionGraph Iterators, which will be described later on in this chapter.
We first show how the programmer can use a runtime system to execute application codes written with these iterators.

In ``amrex/Src/AmrTask/rts_impls``, we implement RTS - a runtime system that can execute asynchronous AMReX applications efficiently on large-scale systems.
RTS is a black box to the application developer as showed in the following code snippet, which is the main function of a typical AMReX application running asynchronously under the control of the runtime system.
The programmer first needs to use the namespace ``perilla``, which covers all the C++ classes for the runtime system.
To execute an AMR program (i.e. object of the Amr class), the programmer can simply create an object of RTS and pass the program object into the ``Iterate`` method. 
The runtime system will iteratively execute coarse time steps until the program completes.
By default RTS links to MPI and Pthreads libraries.
The programmer can also switch to other backends such as UPCXX (1-sided communication model compared to the common 2-sided model in MPI) without changing the application source code.

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

In a few functions of the Amr class, the runtime exposes multiple threads per process.
As a result, the programmer needs to place sufficient memory protection for shared data within the process, e.g. when updating the state data. This multithreaded interface adds some programming cost, but is necessary for mitigating the task scheduling overhead.
To avoid these programming details, the programmer can use built-in iterators, such as fillpatch iterator and task graph iterator that we next discuss.
The API of these iterators is very simple, and the asynchronous code is very similar to the original code using the synchronous multifab iterator (MFIter) described earlier in chapter Basics.


.. toctree::
   :maxdepth: 1

   AsyncIter



