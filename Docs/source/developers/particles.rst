.. _developers-particles:

Particles
=========

Particle containers
-------------------

Particle structures and functions are defined in ``Source/Particles/``. WarpX uses the ``Particle`` class from AMReX for single particles. An ensemble of particles (e.g., a plasma species, or laser particles) is stored as a ``WarpXParticleContainer`` (see description below) in a per-box (and even per-tile on CPU) basis.

.. doxygenclass:: WarpXParticleContainer

Physical species are stored in ``PhysicalParticleContainer``, that derives from ``WarpXParticleContainer``. In particular, the main function to advance all particles in a physical species is ``PhysicalParticleContainer::Evolve`` (see below).

.. doxygenfunction:: PhysicalParticleContainer::Evolve
   :outline:

Finally, all particle species (physical plasma species ``PhysicalParticleContainer``, photon species ``PhotonParticleContainer`` or non-physical species ``LaserParticleContainer``) are stored in ``MultiParticleContainer``. The class ``WarpX`` holds one instance of ``MultiParticleContainer`` as a member variable, called ``WarpX::mypc`` (where `mypc` stands for "my particle containers"):

.. doxygenclass:: MultiParticleContainer

Loop over particles
-------------------

A typical loop over particles reads:

.. code-block:: cpp

  // pc is a std::unique_ptr<WarpXParticleContainer>
  // Loop over MR levels
  for (int lev = 0; lev <= finest_level; ++lev) {
      // Loop over particles, box by box
      for (WarpXParIter pti(*this, lev); pti.isValid(); ++pti) {
          // Do something on particles
          // [MY INNER LOOP]
      }
  }

The innermost step ``[MY INNER LOOP]`` typically calls ``amrex::ParallelFor`` to perform operations on all particles in a portable way. For this reasons, the particle data needs to be converted in plain-old-data structures. The innermost loop in the code snippet above could look like:

.. code-block:: cpp

  // Get Array-Of-Struct particle data, also called data
  // (x, y, z, id, cpu)
  const auto& particles = pti.GetArrayOfStructs();
  // Get Struct-Of-Array particle data, also called attribs
  // (ux, uy, uz, w, Exp, Ey, Ez, Bx, By, Bz)
  auto& attribs = pti.GetAttribs();
  auto& Exp = attribs[PIdx::Ex];
  // [...]
  // Number of particles in this box
  const long np = pti.numParticles();

Link fields and particles?
--------------------------

In WarpX, the loop over boxes through a ``MultiFab`` iterator ``MFIter`` and the loop over boxes through a ``ParticleContainer`` iterator ``WarpXParIter`` are consistent.

On a loop over boxes in a ``MultiFab`` (``MFIter``), it can be useful to access particle data on a GPU-friendly way. This can be done by:

.. code-block:: cpp

  // Index of grid (= box)
  const int grid_id = mfi.index();
  // Index of tile within the grid
  const int tile_id = mfi.LocalTileIndex();
  // Get GPU-friendly arrays of particle data
  auto& ptile = GetParticles(lev)[std::make_pair(grid_id,tile_id)];
  ParticleType* pp = particle_tile.GetArrayOfStructs()().data();
  // Only need attribs (i.e., SoA data)
  auto& soa = ptile.GetStructOfArrays();
  // As an example, let's get the ux momentum
  const ParticleReal * const AMREX_RESTRICT ux = soa.GetRealData(PIdx::ux).data();

On a loop over particles it can be useful to access the fields on the box we are looping over (typically when we use both field and particle data on the same box, for field gather or current deposition for instance). This is done for instance by adding this snippet in ``[MY INNER LOOP]``:

.. code-block:: cpp

  // E is a reference to, say, WarpX::Efield_aux
  // Get the Ex field on the grid
  const FArrayBox& exfab = (*E[lev][0])[pti];
  // Let's be generous and also get the underlying box (i.e., index info)
  const Box& box = pti.validbox();

Main functions
--------------

.. doxygenfunction:: PhysicalParticleContainer::FieldGather

.. doxygenfunction:: PhysicalParticleContainer::PushPX

.. doxygenfunction:: WarpXParticleContainer::DepositCurrent

.. note::
   The current deposition is used both by ``PhysicalParticleContainer`` and ``LaserParticleContainer``, so it is in the parent class ``WarpXParticleContainer``.

Buffers
-------

To reduce numerical artifacts at the boundary of a mesh-refinement patch, WarpX has an option to use buffers: When particles evolve on the fine level, they gather from the coarse level (e.g., ``Efield_cax``, a copy of the ``aux`` data from the level below) if they are located on the fine level but fewer than ``WarpX::n_field_gather_buffer`` cells away from the coarse-patch boundary. Similarly, when particles evolve on the fine level, they deposit on the coarse level (e.g., ``Efield_cp``) if they are located on the fine level but fewer than ``WarpX::n_current_deposition_buffer`` cells away from the coarse-patch boundary.

``WarpX::gather_buffer_masks`` and ``WarpX::current_buffer_masks`` contain masks indicating if a cell is in the interior of the fine-resolution patch or in the buffers. Then, particles depending on this mask in

.. doxygenfunction:: PhysicalParticleContainer::PartitionParticlesInBuffers

.. note::

   Buffers are complex!
