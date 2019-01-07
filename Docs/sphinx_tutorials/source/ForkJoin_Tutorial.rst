.. role:: cpp(code)
   :language: c++

.. role:: fortran(code)
   :language: fortran

******************
Tutorials/Forkjoin
******************

There are two examples in the Tutorials/ForkJoin directory.

``Simple`` demonstrates how to construct a ``ForkJoin`` object, register
some MultiFabs, and invoke the fork_join operation.

``MLMG`` demonstrates how to do more advanced fork-join operations,
including nested fork-joins, heterogeneous tasks, reusing ``ForkJoin``
objects, and customized MultiFab component splitting.

**General Concepts**
====================

An AMReX program consists of a set of MPI ranks cooperating together on
distributed data.
Typically, all of the ranks in a job compute in a bulk-synchronous,
data-parallel fashion, where every rank does the same sequence of
operations, each on different parts of the distributed data.

The ForkJoin functionality described here allows the user to divide the
job's MPI ranks into subgroups (i.e. `fork`) and assign each subgroup
an independent task to compute in parallel with each other.
After all of the forked child tasks complete, they synchronize
(i.e. `join`), and the parent task continues execution as before.

.. figure:: fork_join_tasks.png
   :scale: 50 %
   :alt: Fork-Join Tasks

   Example of a fork-join operation where the parent task's MPI ranks are
   split into two independent child tasks that execute in parallel and
   then join to resume execution of the parent task.

The fork-join operation can also be invoked in a nested fashion,
creating a hierarchy of fork-join operations, where each fork further
subdivides the ranks of a task into child tasks.
This approach enables heterogeneous computation and reduces the strong
scaling penalty for operations with less inherent parallelism or with
large communication overheads.

.. figure:: nested_fork_join_tasks.png
   :scale: 50 %
   :alt: Nested Fork-Join Tasks

   Example of nested fork-join operations where a child task is further
   split into more subtasks.

The fork-join operation is accomplished by a) redistributing MultiFab
data so that all of the data in each (registered) MultiFab is visible
to ranks within a subtask, and b) dividing the root MPI communicator
into sub-communicators so that each subgroup of ranks in a tasks will
only synchronize with each other during subtask collectives (e.g. for
``MPI_Allreduce``).
When the program starts, all of the ranks in the MPI communicator are
in the root task.

**ForkJoin/Simple**
===================

The main function in this tutorial is in ``MyTest.cpp:runTest()``.
It does the following things:

1. Create a ``ForkJoin`` object: the constructor takes the number of
tasks to split the calling (in this case, root) task.
This version of the constructor will divide the ranks in the calling
(parent) task evenly across the spawned (child) tasks.
To allow uneven distribution of ranks across tasks, there are other
versions of the ``ForkJoin`` constructor that allow the user to specify
the number (or percent) of ranks to include in each of the subtasks.

2. Set the verbosity flag and task output directory:
There are member functions of the ``ForkJoin`` object to set each of
these behaviors.
Both of these calls may be omitted to accept default values.

3. Register three MultiFab data structures:
The ``ForkJoin`` object needs to know what data will be utilized within the
spawned subtasks and how they will be accessed.
For each MultiFab that will be accessed within the subtasks, there are
two main parameters that need to be specified: Strategy and Intent.
``Strategy`` describes whether the MultiFab will be ``duplicate`` across
all tasks, ``split`` (component-wise) across the subtasks, or accessed
in only a ``single`` subtask.
``Intent`` describes whether the data is an input and/or output to the
forked subtasks, and controls whether the data is copied in and/or out
of the subtask from the calling task.

.. figure:: mf_remap_hires.png
   :scale: 50 %
   :alt: Examples of how to register MultiFabs

   Examples of how a MultiFab can be registered for a fork-join operation
   with varying Strategy and Intent.

4. Invoke the fork_join operation by calling ``myFunction`` in every task:
The fork_join function launches the passed function (or lambda) on
all of the spawned tasks.
The passed function must take a single argument: a reference to the
managing ``ForkJoin`` object, which can be queried for the subtask's ID,
references to the registered MultiFabs, and other metadata such as the
component bounds of a registered MultiFab.
The tutorial's ``myFunction`` demonstrates these capabilities.

**ForkJoin/MLMG**
=================


