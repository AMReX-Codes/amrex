.. role:: cpp(code)
   :language: c++

.. role:: fortran(code)
   :language: fortran

Amr Class
=========

The :cpp:`Amr` class is designed to manage parts of the computation which do
not belong on a single level, like establishing and updating the hierarchy of
levels, global timestepping, and managing the different :cpp:`AmrLevel`\ s.
Most likely you will not need to derive any classes from :cpp:`Amr`. Our
mature application codes use this base class without any derived classes.

One of the most important data members is an array of :cpp:`AmrLevel`\ s - the
:cpp:`Amr` class calls many functions from the :cpp:`AmrLevel` class to do
things like advance the solution on a level, compute a time step to be used for
a level, etc.

AmrLevel Class
==============

Pure virtual functions include:

   -  :cpp:`computeInitialDt` Compute an array of time steps for each level of
      refinement. Called at the beginning of the simulation.

   -  :cpp:`computeNewDt` Compute an array of time steps for each level of
      refinement. Called at the end of a coarse level advance.

   -  :cpp:`advance` Advance the grids at a level.

   -  :cpp:`post_timestep` Work after a time step at a given level. The
      AMR synchronization is usually done here.

   -  :cpp:`post_regrid` Work after regridding, e.g., redistributing
      particles.

   -  :cpp:`post_init` Work after initialization, e.g., AMR synchronization.

   -  :cpp:`initData` Initialize the data on a given level at the beginning of
      the simulation.

   -  :cpp:`init` There are two versions of this function used to initialize
      data on a level during regridding. One version is specifically for the
      case where the level did not previously exist (a newly created refined
      level).

   -  :cpp:`errorEst` Perform the tagging at a level for refinement.

StateData
---------

The most important data managed by the :cpp:`AmrLevel` is an array of
:cpp:`StateData`, which holds the scalar fields, etc., in the boxes that
together make up the level.

:cpp:`StateData` is a class that essentially holds a pair of MultiFabs: one at
the old time and one at the new time. AMReX knows how to interpolate in time
between these states to get data at any intermediate point in time. The main
data that we care about in our application codes (such as the fluid state)
will be stored as :cpp:`StateData`.  Essentially, data is made :cpp:`StateData`
if we need it to be stored in checkpoints/plotfiles, and/or we want it to be
automatically interpolated when we refine.  An :cpp:`AmrLevel` stores an array
of :cpp:`StateData` (in a C++ array called :cpp:`state`). We index this array
using integer keys (defined via an :cpp:`enum` in, e.g., ``AmrLevelWave.H`` in
the WaveAmrLevel tutorial):

.. highlight:: c++

::

    enum StateType { State_Type = 0, NUM_STATE_TYPE };

In the tutorial code, we use the function :cpp:`AmrLevelWave::variableSetUp` to
tell our simulation about the :cpp:`StateData` (e.g., how many variables, ghost
cells, nodality, etc.). Note that if you have more than one :cpp:`StateType`,
each of the different :cpp:`StateData` carried in the state array can have
different numbers of components, ghost cells, boundary conditions, etc.
This is the main reason we separate all this data into separate :cpp:`StateData`
objects collected together in an indexable array.

LevelBld Class
==============

The :cpp:`LevelBld` class is a pure virtual class for defining variable types
and attributes. To more easily understand its usage, refer to the derived
class, :cpp:`LevelBldWave` in the WaveAmrLevel tutorial. The :cpp:`variableSetUp` and
:cpp:`variableCleanUp` functions are implemented, and in this tutorial they call routines in
the :cpp:`AmrLevelWave` class, e.g.,

.. highlight:: c++

::

    void
    AmrLevelWave::variableSetUp ()
    {
        read_params();

        // One cell-centered StateData with two components and the ghost cells
        // of the stencil. Fine-level ghost cells at the coarse/fine boundary
        // are interpolated from the coarse level with a fourth-order
        // interpolater, to match the fourth-order spatial discretization.
        desc_lst.addDescriptor(State_Type, IndexType::TheCellType(),
                               StateDescriptor::Point, nghost, ncomp,
                               &cell_quartic_interp);

        int lo_bc[AMREX_SPACEDIM] = {AMREX_D_DECL(BCType::ext_dir,    // external Dirichlet
                                                  BCType::int_dir,    // periodic
                                                  BCType::int_dir)};  // periodic
        int hi_bc[AMREX_SPACEDIM] = {AMREX_D_DECL(BCType::ext_dir,
                                                  BCType::int_dir,
                                                  BCType::int_dir)};
        Vector<BCRec> bcs(ncomp, BCRec(lo_bc, hi_bc));

        StateDescriptor::BndryFunc bndryfunc(wave_bcfill);
        bndryfunc.setRunOnGPU(true);

        desc_lst.setComponent(State_Type, 0, {"u", "v"}, bcs, bndryfunc);
    }

We see how to define the :cpp:`StateType`, including nodality, whether or not
we want the variable to represent a point in time or an interval over time
(useful for returning the time associated with data), the number of ghost
cells, number of components, and the interlevel interpolation (See
AMReX_Interpolator for various interpolation types). We also see how to specify
physical boundary functions by providing a function (in this case,
:cpp:`wave_bcfill`, defined in ``AmrLevelWave.cpp`` in the tutorial source
code, which fills the external Dirichlet boundaries in the :math:`x`
direction).

Examples
========

The `WaveAmrLevel`_ tutorial is a small code built on :cpp:`Amr` and
:cpp:`AmrLevel`. It contains the :cpp:`AmrLevelWave` class, derived from the
base :cpp:`AmrLevel` class, and the :cpp:`LevelBldWave` class, derived from
the base :cpp:`LevelBld` class as described above. It advances the state
with :cpp:`AmrLevel::RK`, which provides Runge-Kutta time stepping with
subcycling.

The `CNS`_ tutorial is a complete compressible flow code built on the same
classes, with refluxing, derived variables for plotfiles and GPU support.
It is meant to be copied as the starting point of an application.

.. _`WaveAmrLevel`: https://amrex-codes.github.io/amrex/tutorials_html/AMR.html#amr-waveamrlevel

.. _`CNS`: https://amrex-codes.github.io/amrex/tutorials_html/Applications.html#applications-cns

The pseudocode for the main program of such a code is given below.

.. highlight:: c++

::

    /* AmrLevel Pseudocode */
    main()
      Amr amr(&level_bld);
      amr.init()
      loop {
        amr.coarseTimeStep()
          /* compute dt */
          timeStep()
            amr_level[level]->advance()
            /* call timeStep r times for next-finer level */
            amr_level[level]->post_timestep() // AMR synchronization
          postCoarseTimeStep()
          /* write plotfile and checkpoint */
      }
      /* write final plotfile and checkpoint */

Particles
=========

``amrex/Tests/Amr/Advection_AmrLevel`` is another :cpp:`AmrLevel` code, with an
option to turn on passively advected particles. Its
``Exec/SingleVortex/GNUmakefile`` builds with ``USE_PARTICLES = TRUE``, and
the ``inputs.tracers`` file there turns the particles on with
``adv.do_tracers = 1``.
When you run the code, within each plotfile directory there will be a
subdirectory called "Tracer". These can be visualized using either yt or
ParaView (refer to the chapter on :ref:`Chap:Visualization`).
