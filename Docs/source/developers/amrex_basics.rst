.. _developers-amrex-basics:

AMReX basics (excessively basic)
================================

WarpX is built on the Adaptive Mesh Refinement (AMR) library `AMReX <https://github.com/AMReX-Codes/amrex>`__. This section provides a very sporadic description of the main AMReX classes and concepts relevant for WarpX, that can serve as a reminder. Please read the AMReX basics `doc page <https://amrex-codes.github.io/amrex/docs_html/Basics.html>`__, of which this section is largely inspired.

* ``amrex::Box``: Dimension-dependent lower and upper indices defining a rectangular volume in 3D (or surface in 2D) in the index space. ``Box`` is a lightweight meta-data class, with useful member functions.

* ``amrex::BoxArray``: Collection of ``Box`` on a single AMR level. The information of which MPI rank owns which ``Box`` in a ``BoxArray`` is in ``DistributionMapping``.

* ``amrex::FArrayBox``: Fortran-ordered array of floating-point ``amrex::Real`` elements defined on a ``Box``. A ``FArrayBox`` can represent scalar data or vector data, with ``ncomp`` components.

* ``amrex::MultiFab``: Collection of `FAB` (= ``FArrayBox``) on a single AMR level, distributed over MPI ranks. The concept of `ghost cells` is defined at the ``MultiFab`` level.

* ``amrex::ParticleContainer``: A collection of particles, typically for particles of a physical species. Particles in a ``ParticleContainer`` are organized per ``Box``. Particles in a ``Box`` are organized per tile (this feature is off when running on GPU). Particles within a tile are stored in several structures, each being contiguous in memory: (i) an Array-Of-Struct (AoS) (often called `data`, they are the 3D position, the particle ID and the index of the CPU owning the particle), where the Struct is an ``amrex::Particle`` and (ii) Struct-Of-Arrays (SoA) for extra variables (often called ``attribs``, in WarpX they are the momentum, field on particle etc.).

The simulation domain is decomposed in several ``Box``, and each MPI rank owns (and performs operations on) the fields and particles defined on a few of these ``Box``, but has the metadata of all of them. For convenience, AMReX provides iterators, to easily iterate over all ``FArrayBox`` (or even tile-by-tile, optionally) in a ``MultiFab`` own by the MPI rank (``MFIter``), or over all particles in a ``ParticleContainer`` on a per-box basis (``ParIter``, or its derived class ``WarpXParIter``). These are respectively done in loops like:

.. code-block:: cpp

   // mf is a pointer to MultiFab
   for ( amrex::MFIter mfi(mf, false); mfi.isValid(); ++mfi ) { ... }

and

.. code-block:: cpp

   // *this is a pointer to a ParticleContainer
   for (WarpXParIter pti(*this, lev); pti.isValid(); ++pti) { ... }

When looping over ``FArrayBox`` in a ``MultiFab``, the iterator provides functions to retrieve the metadata of the ``Box`` on which the ``FAB`` is defined (``MFIter::box()``, ``MFIter::tilebox()`` or variations) or the particles defined on this ``Box`` (``ParIter::GetParticles()``).
