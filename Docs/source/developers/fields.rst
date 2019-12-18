Fields
======

The main fields are the electric field ``Efield``, the magnetic field ``Bfield``, the current density ``current`` and the charge density ``rho``. When a divergence-cleaner is used, another field ``F`` (containing ``div(Efield)-rho``).

Due the AMR strategy used in WarpX (see section :doc:`../../theory/amr` for a complete description), each field on a given refinement level ``lev`` (except for the coarsest `0`) is defined on:

* **the fine patch** (suffix ``_fp``, the actual resolution on ``lev``).

* **the coarse patch** (suffix ``_cp``, same physical domain with the resolution of MR level ``lev-1``).

* **the auxiliary grid** (suffix ``_aux``, same resolution as ``_fp``), from which the fields are gathered from the grids to particle positions. For this reason. only `E` and `B` are defined on this ``_aux`` grid (not the current density or charge density).

* In some conditions, i.e., when buffers are used for the field gather (for numerical reasons), a **copy of E and B on the auxiliary grid** ``_aux`` **of the  level below** ``lev-1`` is stored in fields with suffix ``_cax`` (for `coarse aux`).

As an example, the structures for the electric field are ``Efield_fp``, ``Efield_cp``, ``Efield_aux`` (and optionally ``Efield_cax``).

Declaration
-----------

All the fields described above are public members of class ``WarpX``, defined in ``WarpX.H``. They are defined as an ``amrex::Vector`` (over MR levels) of ``std::array`` (for the 3 spatial components `Ex`, `Ey`, `Ez`) of ``std::unique_ptr`` of ``amrex::MultiFab``, i.e.,
::

  amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3 > > Efield_fp;

Hence, `Ex` on MR level ``lev`` is a pointer to an ``amrex::MultiFab``. The other fields are organized in the same way.

Allocation and initialization
-----------------------------

The ``MultiFab`` constructor (for, e.g., `Ex` on level `lev`) is called in ``WarpX::AllocLevelMFs``.

By default, the ``MultiFab`` are set to ``0`` at initialization. They can be assigned a different value in ``WarpX::InitLevelData``.

.. note::

   Add info on staggering, guard cells and domain decomposition. Synchronize with section ``initialization``.

Field solver
------------

The field solver is performed in ``WarpX::EvolveE`` for the electric field and ``WarpX::EvolveB`` for the magnetic field, called from ``WarpX::OneStep_nosub`` in ``WarpX::EvolveEM``. This section describes the FDTD field push (the PSATD field push should be added later). It is implemented in ``Source/FieldSolver/WarpXPushFieldsEM.cpp``.

As all cell-wise operation, the field push is done as follows (this is split in multiple functions in the actual implementation to aboid code duplication)::

  // Loop over MR levels
  for (int lev = 0; lev <= finest_level; ++lev) {
      // Get pointer to MultiFab Ex on level lev
      MultiFab* Ex = Efield_fp[lev][0].get();
      // Loop over boxes (or tiles if not on GPU)
      for ( MFIter mfi(*Ex, TilingIfNotGPU()); mfi.isValid(); ++mfi ) {
          // Apply field solver on the FAB
      }
  }

The innermost step ``// Apply field solver on the FAB`` could be done with 3 nested ``for`` loops for the 3 dimensions (in 3D). However, for portability reasons (see section :doc:`./portability`), this is done in two steps: (i) extract AMReX data structures into plain-old-data simple structures, and (ii) call a general ``ParallelFor`` function (translated into nested loops on CPU or a kernel launch on GPU, for instance)::

  // Get Box corresponding to the current MFIter
  const Box& tex  = mfi.tilebox(Ex_nodal_flag);
  // Extract the FArrayBox into a simple structure, for portability
  Array4<Real> const& Exfab = Ex->array(mfi);
  // Loop over cells and perform stencil operation
  amrex::ParallelFor(tex,
      [=] AMREX_GPU_DEVICE (int j, int k, int l)
      {
          warpx_push_ex_yee(...);
      }
  );

Function ``warpx_push_ex_yee`` performs the FDTD stencil operation on a single cell. It is implemented in ``Source/FieldSolver/WarpX_K.H`` (where ``_K`` stands for kernel).

Guard cells exchanges
---------------------

Communications are mostly handled in ``Source/Parallelization/``.

For E and B guard cell **exchanges**, the main functions are variants of ``amrex::FillBoundary(amrex::MultiFab, ...)`` (or ``amrex::MultiFab::FillBoundary(...)``) that fill guard cells of all ``amrex::FArrayBox`` in an ``amrex::MultiFab`` with valid cells of corresponding ``amrex::FArrayBox`` neighbors of the same ``amrex::MultiFab``. There are a number of ``FillBoundaryE``, ``FillBoundaryB`` etc. Under the hood, ``amrex::FillBoundary`` calls ``amrex::ParallelCopy``, which is also sometimes directly called in WarpX. Most calls a

For the current density, the valid cells of neighboring ``MultiFabs`` are accumulated (added) rather than just copied. This is done using ``amrex::MultiFab::SumBoundary``, and mostly located in ``Source/Parallelization/WarpXSumGuardCells.H``.

Interpolations for MR
---------------------

This is mostly implemented in ``Source/Parallelization``, see the following functions (you may complain to the authors if the documentation is empty)

.. doxygenfunction:: WarpX::SyncCurrent

.. doxygenfunction:: WarpX::interpolateCurrentFineToCoarse

.. doxygenfunction:: WarpX::RestrictCurrentFromFineToCoarsePatch

.. doxygenfunction:: WarpX::AddCurrentFromFineLevelandSumBoundary

Filter
------

General functions for filtering can be found in ``Source/Filter/``, where the main ``Filter`` class is defined (see below). All filters (so far there are two of them) in WarpX derive from this class.

.. doxygenclass:: Filter

Bilinear filter
~~~~~~~~~~~~~~~

See

Godfrey's anti-NCI filter for FDTD simulations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


Buffers
-------

.. note::

   Section empty!
