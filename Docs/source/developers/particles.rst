Particles
=========

Particle structures and functions are defined in ``Source/Particles/``. WarpX uses the ``Particle`` class from AMReX for single particles. The main class for ensembles of particles is ``WarpXParticleContainer``.

.. doxygenclass:: WarpXParticleContainer

``PhysicalParticleContainer``, ``PhotonParticleContainer`` and ``LaserParticleContainer`` inherit from ``WarpXParticleContainer``. ``RigidInjectedParticleContainer`` inherits from ``PhysicalParticleContainer``. The main function ``WarpXParticleContainer::Evolve`` is overridden by derived classes (see below).

.. doxygenfunction:: PhysicalParticleContainer::Evolve
