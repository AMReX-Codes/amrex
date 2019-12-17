Particles
=========

Particle structures and functions are defined in ``Source/Particles/``. WarpX uses the ``Particle`` class from AMReX for single particles. An ensemble of particles (e.g., a plasma species, or laser particles) is stored as a ``WarpXParticleContainer`` (see description below) in a per-box (and even per-tile on CPU) basis.

.. doxygenclass:: WarpXParticleContainer

Physical species are stored in ``PhysicalParticleContainer``, that derives from ``WarpXParticleContainer``. In particular, the main function to advance all particles in a physical species is ``PhysicalParticleContainer::Evolve`` (see below).

.. doxygenfunction:: PhysicalParticleContainer::Evolve
   :outline:

Finally, all particle species (physical plasma species ``PhysicalParticleContainer``, photon species ``PhotonParticleContainer`` or non-physical species ``LaserParticleContainer``) are stored in ``MultiParticleContainer``. The class ``WarpX`` holds one instance of ``MultiParticleContainer`` as a member variable, called ``WarpX::mypc`` (where `mypc` stands for "my particle containers"):

.. doxygenclass:: MultiParticleContainer

