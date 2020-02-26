.. _developers-fields:

Fields
======

.. note::

   Add info on staggering and domain decomposition. Synchronize with section ``initialization``.

The main fields are the electric field ``Efield``, the magnetic field ``Bfield``, the current density ``current`` and the charge density ``rho``. When a divergence-cleaner is used, we add another field ``F`` (containing :math:`\vec \nabla \cdot \vec E - \rho`).

Due the AMR strategy used in WarpX (see section :ref:`Theory: AMR <theory-amr>` for a complete description), each field on a given refinement level ``lev`` (except for the coarsest ``0``) is defined on:

* **the fine patch** (suffix ``_fp``, the actual resolution on ``lev``).

* **the coarse patch** (suffix ``_cp``, same physical domain with the resolution of MR level ``lev-1``).

* **the auxiliary grid** (suffix ``_aux``, same resolution as ``_fp``), from which the fields are gathered from the grids to particle positions. For this reason. only ``E`` and ``B`` are defined on this ``_aux`` grid (not the current density or charge density).

* In some conditions, i.e., when buffers are used for the field gather (for numerical reasons), a **copy of E and B on the auxiliary grid** ``_aux`` **of the  level below** ``lev-1`` is stored in fields with suffix ``_cax`` (for `coarse aux`).

As an example, the structures for the electric field are ``Efield_fp``, ``Efield_cp``, ``Efield_aux`` (and optionally ``Efield_cax``).

Declaration
-----------

All the fields described above are public members of class ``WarpX``, defined in ``WarpX.H``. They are defined as an ``amrex::Vector`` (over MR levels) of ``std::array`` (for the 3 spatial components :math:`E_x`, :math:`E_y`, :math:`E_z`) of ``std::unique_ptr`` of ``amrex::MultiFab``, i.e.:

  amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3 > > Efield_fp;

Hence, ``Ex`` on MR level ``lev`` is a pointer to an ``amrex::MultiFab``. The other fields are organized in the same way.

Allocation and initialization
-----------------------------

The ``MultiFab`` constructor (for, e.g., ``Ex`` on level ``lev``) is called in ``WarpX::AllocLevelMFs``.

By default, the ``MultiFab`` are set to ``0`` at initialization. They can be assigned a different value in ``WarpX::InitLevelData``.

Field solver
------------

The field solver is performed in ``WarpX::EvolveE`` for the electric field and ``WarpX::EvolveB`` for the magnetic field, called from ``WarpX::OneStep_nosub`` in ``WarpX::EvolveEM``. This section describes the FDTD field push. It is implemented in ``Source/FieldSolver/FiniteDifferenceSolver/``.

As all cell-wise operation, the field push is done as follows (this is split in multiple functions in the actual implementation to aboid code duplication)
:

.. code-block:: cpp

   // Loop over MR levels
   for (int lev = 0; lev <= finest_level; ++lev) {
      // Get pointer to MultiFab Ex on level lev
      MultiFab* Ex = Efield_fp[lev][0].get();
      // Loop over boxes (or tiles if not on GPU)
      for ( MFIter mfi(*Ex, TilingIfNotGPU()); mfi.isValid(); ++mfi ) {
          // Apply field solver on the FAB
      }
  }

The innermost step ``// Apply field solver on the FAB`` could be done with 3 nested ``for`` loops for the 3 dimensions (in 3D). However, for portability reasons (see section :ref:`Developers: Portability <developers-portability>`), this is done in two steps: (i) extract AMReX data structures into plain-old-data simple structures, and (ii) call a general ``ParallelFor`` function (translated into nested loops on CPU or a kernel launch on GPU, for instance):

.. code-block:: cpp

  // Get Box corresponding to the current MFIter
  const Box& tex  = mfi.tilebox(Ex_nodal_flag);
  // Extract the FArrayBox into a simple structure, for portability
  Array4<Real> const& Exfab = Ex->array(mfi);
  // Loop over cells and perform stencil operation
  amrex::ParallelFor(tex,
      [=] AMREX_GPU_DEVICE (int j, int k, int l)
      {
          Ex(i, j, k) += c2 * dt * (
              - T_Algo::DownwardDz(By, coefs_z, n_coefs_z, i, j, k)
              + T_Algo::DownwardDy(Bz, coefs_y, n_coefs_y, i, j, k)
              - PhysConst::mu0 * jx(i, j, k) );
      }
  );

where ``T_Algo::DownwardDz`` and ``T_Algo::DownwardDy`` represent the discretized derivative
for a given algorithm (represented by the template parameter ``T_Algo``). The available
discretization algorithms can be found in ``Source/FieldSolver/FiniteDifferenceSolver/FiniteDifferenceAlgorithms``.

Guard cells exchanges
---------------------

Communications are mostly handled in ``Source/Parallelization/``.

For E and B guard cell **exchanges**, the main functions are variants of ``amrex::FillBoundary(amrex::MultiFab, ...)`` (or ``amrex::MultiFab::FillBoundary(...)``) that fill guard cells of all ``amrex::FArrayBox`` in an ``amrex::MultiFab`` with valid cells of corresponding ``amrex::FArrayBox`` neighbors of the same ``amrex::MultiFab``. There are a number of ``FillBoundaryE``, ``FillBoundaryB`` etc. Under the hood, ``amrex::FillBoundary`` calls ``amrex::ParallelCopy``, which is also sometimes directly called in WarpX. Most calls a

For the current density, the valid cells of neighboring ``MultiFabs`` are accumulated (added) rather than just copied. This is done using ``amrex::MultiFab::SumBoundary``, and mostly located in ``Source/Parallelization/WarpXSumGuardCells.H``.

Interpolations for MR
---------------------

This is mostly implemented in ``Source/Parallelization``, see the following functions (you may complain to the authors if the documentation is empty)

.. doxygenfunction:: WarpX::SyncCurrent

.. doxygenfunction:: interpolateCurrentFineToCoarse

.. doxygenfunction:: WarpX::RestrictCurrentFromFineToCoarsePatch

.. doxygenfunction:: WarpX::AddCurrentFromFineLevelandSumBoundary

Filter
------

General functions for filtering can be found in ``Source/Filter/``, where the main ``Filter`` class is defined (see below). All filters (so far there are two of them) in WarpX derive from this class.

.. doxygenclass:: Filter

Bilinear filter
~~~~~~~~~~~~~~~

The multi-pass bilinear filter (applied on the current density) is implemented in ``Source/Filter/``, and class ``WarpX`` holds an instance of this class in member variable ``WarpX::bilinear_filter``. For performance reasons (to avoid creating too many guard cells), this filter is directly applied in communication routines, see

.. doxygenfunction:: WarpX::AddCurrentFromFineLevelandSumBoundary

and

.. doxygenfunction:: WarpX::ApplyFilterandSumBoundaryJ

Godfrey's anti-NCI filter for FDTD simulations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This filter is applied on the electric and magnetic field (on the auxiliary grid) to suppress the Numerical Cherenkov Instability when running FDTD. It is implemented in ``Source/Filter/``, and there are two different stencils, one for ``Ex``, ``Ey`` and ``Bz`` and the other for ``Ez``, ``Bx`` and ``By``.

.. doxygenclass:: NCIGodfreyFilter

The class ``WarpX`` holds two corresponding instances of this class in member variables ``WarpX::nci_godfrey_filter_exeybz`` and ``WarpX::nci_godfrey_filter_bxbyez``. It is a 9-point stencil (is the ``z`` direction only), for which the coefficients are computed using tabulated values (depending on dz/dx) in ``Source/Utils/NCIGodfreyTables.H``, see variable ``table_nci_godfrey_galerkin_Ex_Ey_Bz``. The filter is applied in ``PhysicalParticleContainer::Evolve``, right after field gather and before particle push, see

.. doxygenfunction:: PhysicalParticleContainer::applyNCIFilter
