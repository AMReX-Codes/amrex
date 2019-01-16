.. role:: cpp(code)
   :language: c++

.. role:: fortran(code)
   :language: fortran

Tutorials/GPU
==========================

The tutorials in ``amrex/Tutorials/GPU`` demonstrate the implementation
of AMReX's GPU toolkit as well as provide GPU ported versions of CPU
tutorials to help applications convert to GPUs. 

**Your first AMReX GPU application**
------------------------------------

This is a step-by-step guide to preparing, compiling and running your first
AMReX GPU program.  This guide will use ``Tutorials/GPU/Launch``,
and instructions will focus on ORNL's Summit system:

1. Before compiling, the ``pgi`` and ``cuda`` software must be available.  On
   ORNL systems, this is done with ``module load pgi cuda``.  However, Cuda
   versions 9.2.x are not compatible with AMReX.  So, it is important to check
   whether your system is using 9.2.x as its default. To check, type 
   ``module list``.  If Cuda is the wrong version, use ``module avail cuda`` 
   to find a compatible version and swap using 
   ``module swap cuda [good module name]``.

2. Go to ``Tutorials/GPU/Launch`` to compile the executable.  Compile with
   ``make USE_CUDA=TRUE COMP=pgi USE_MPI=TRUE USE_OMP=FALSE``. This should
   result in an executable: ``main3d.pgi.MPI.CUDA.ex``.  

3. On Summit systems, this executable can be submitted by using one of the run
   scripts found in ``Tutorials/GPU``.  ``run.script`` can be used to run on
   Summitdev, and ``run.summit`` can be used for Summit.  To change the number
   of ranks and GPUs used in the simulation, change the number of resource sets,
   ``n`` in the ``jsrun`` line.  Set the inputs line to an empty string for this
   case. 

When ready, submit the job (on Summit: ``bsub run.script``).  Congratulations!
You've accelerated AMReX using GPUs! 

**Launch**
----------

Launch shows multiple examples of how GPU work can be offloaded using the tools
available in AMReX. It includes examples of multiple AMReX macro launch methods,
launching a Fortran function using CUDA and launching work using OpenACC and 
OpenMP offloading. This tutorial will be regularly updated with AMReX's 
preferred GPU launch methodologies.

**HeatEquation_EX1_C**
----------------------

HeatEquation is a direct GPU port of the ``Tutorials/Basic/HeatEquation_EX1_C``
tutorial that solves the 2D or 3D heat equation on a domain-decomposed mesh. It
offloads the phi :cpp:`Multifab` initialization, flux computation and phi update
kernels to the GPU using AMReX CUDA launch macros and has converted all Fortran
loops into C++ inlined functions. 

**CNS**
-------

CNS is a direct GPU port of the ``Tutorials/EB/CNS`` tutorial.

**AmrCore**
-----------

AmrCore is a direct GPU port of the ``Tutorials/Amr/Advection_AmrCore`` tutorial
that advects a single scalar field with a velocity field specified on faces, using
strategies similar to HeatEquation and CNS.

