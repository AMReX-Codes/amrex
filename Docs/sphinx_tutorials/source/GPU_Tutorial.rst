.. role:: cpp(code)
   :language: c++

.. role:: fortran(code)
   :language: fortran

Tutorials/GPU
==========================

The tutorials in ``amrex/Tutorials/GPU`` demonstrate the implementation
of AMReX's GPU toolkit as well as provide GPU versions of other CPU
tutorials for comparison to help applications convert to GPUs. 

**Your first AMReX GPU application**
------------------------------------

This is a step-by-step guide to preparing, compiling and running your first
AMReX GPU program.  This guide will use ``Tutorials/GPU/HeatEquation_EX1_C``,
and instructions will focus on ORNL's Summit and Summitdev systems.

Before compiling, the ``pgi`` and ``cuda`` software must be available.  On
ORNL systems, this is done with ``module load pgi cuda``.  However, Cuda
versions 9.2.x are not compatible with AMReX.  So, it is important to check
whether your system is using 9.2.x as its default. To check, type 
``module list``.  If Cuda is the wrong version, use ``module avail cuda`` 
to find a compatible version and swap using 
``module swap cuda [good module name]``.

Go to ``Tutorials/GPU/HeatEquation_EX1_C/Exec/SingleVortex`` to compile the
executable.  Compile with ``USE_CUDA=TRUE``, ``COMP=pgi``, ``USE_MPI=TRUE``
and ``USE_OMP=FALSE``.  This should result in an executable: 
``main3d.pgi.MPI.CUDA.ex``.  

On ORNL's systems, this executable can be ran by using one of the run scripts
found in ``Tutorials/GPU``. ``run.script`` can be used to run on Summitdev,
and ``run.summit`` can be used for Summit.  To change the number of ranks and
GPUs used in the simulation, change the number of resource sets, ``n`` in the
``jsrun`` linei. Also, set the inputs line to take the appopriate input file
based on the dimensionality of your build. 

When ready, submit the job (on ORNL: ``bsub run.script``). Congratulations!
You've accelerated AMReX using GPUs! 

**Launch**
----------

**HeatEquation_EX1_C**
----------------------

**CNS**
-------

**AmrCore**
-----------

