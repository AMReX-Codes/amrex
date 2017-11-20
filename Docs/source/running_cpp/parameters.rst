Input parameters
================

.. warning::

   This section is currently in development.


Overall simulation parameters
-----------------------------

* ``max_step`` (`integer`)
    The number of PIC cycles to perform.


* ``warpx.gamma_boost`` (`float`)
    The Lorentz factor of the boosted frame in which the simulation is run.
    (The corresponding Lorentz transformation is assumed to be along ``warpx.boost_direction``.)

    When using this parameter, some of the input parameters are automatically
    converted to the boosted frame. (See the corresponding documentation of each
    input parameters.)

    .. note::

        For now, only the laser parameters will be converted.

* ``warpx.boost_direction`` (`3 floats in 3D`)
    The coordinates of a vector that points in the propagation direction of
    the Lorentz transform. The norm of this vector is unimportant, only its direction matters.

Setting up the field mesh
-------------------------

* ``amr.n_cell`` (`2 integers in 2D`, `3 integers in 3D`)
    The number of grid points along each direction (on the **coarsest level**)

* ``amr.max_level`` (`integer`)
    When using mesh refinement, the number of refinement levels that will be used.

    Use 0 in order to disable mesh refinement.

* ``geometry.is_periodic`` (`2 integers in 2D`, `3 integers in 3D`)
    Whether the boundary conditions are periodic, in each direction.

    For each direction, use 1 for periodic conditions, 0 otherwise.

* ``geometry.prob_lo`` and ``geometry.prob_hi`` (`2 floats in 2D`, `3 integers in 3D`; in meters)
    The extent of the full simulation box. This box is rectangular, and thus its
    extent is given here by the coordinates of the lower corner (``geometry.prob_lo``) and
    upper corner (``geometry.prob_hi``).

* ``warpx.fine_tag_lo`` and ``warpx.fine_tag_hi`` (`2 floats in 2D`, `3 integers in 3D`; in meters)
    **When using static mesh refinement with 1 level**, the extent of the refined patch.
    This patch is rectangular, and thus its extent is given here by the coordinates
    of the lower corner (``warpx.fine_tag_lo``) and upper corner (``warpx.fine_tag_hi``).

Distribution across MPI ranks and parallelization
-------------------------------------------------


* ``amr.max_grid_size`` (`integer`)
    Maximum allowable size of each **subdomain**
    (expressed in number of grid points, in each direction).
    Each subdomain has its own ghost cells, and can be handled by a
    different MPI rank ; several OpenMP threads can work simultaneously on the
    same subdomain.

    If ``max_grid_size`` is such that the total number of subdomains is
    **larger** that the number of MPI ranks used, than some MPI ranks
    will handle several subdomains, thereby providing additional flexibility
    for **load balancing**.

    When using mesh refinement, this number applies to the subdomains
    of the coarsest level, but also to any of the finer level.

* ``warpx.load_balance_int`` (`integer`)
    How often WarpX should try to redistribution the work across MPI ranks,
    in order to have better load balancing (expressed in number of PIC cycles
    inbetween two consecutive attempts at redistributing the work).
    Use 0 to disable load_balancing.

    When performing load balancing, WarpX measures the wall time for
    computational parts of the PIC cycle. It then uses this data to decide
    how to redistribute the subdomains across MPI ranks. (Each subdomain
    is unchanged, but its owner is changed in order to have better performance.)
    This relies on each MPI rank handling several (in fact many) subdomains
    (see ``max_grid_size``).


Particle initialization
-----------------------

Laser initialization
--------------------

* ``warpx.use_laser`` (`0 or 1`)
    Whether to activate the injection of a laser pulse in the simulation

* ``laser.profile`` (`string`)
    The spatio-temporal shape of the laser. The options that are currently
    implemented are:

    - ``"Gaussian"``: The transverse and longitudinal profiles are Gaussian.
    - ``"Harris"``: The transverse profile is Gaussian, but the longitudinal profile is given by the Harris function (see ``laser.profile_duration`` for more details)

* ``laser.e_max`` (`float` ; in V/m)
    Peak amplitude of the laser field.

    For a laser with a wavelength :math:`\lambda = 0.8\,\mu m`, the peak amplitude
    is related to :math:`a_0` by:

    .. math::

        E_{max} = a_0 \frac{2 \pi m_e c}{e\lambda} = a_0 \times (4.0 \cdot 10^{12} \;V.m^{-1})

    When running a **boosted-frame simulation**, provide the value of ``laser.e_max``
    in the laboratory frame, and use ``warpx.gamma_boost`` to automatically
    perform the conversion to the boosted frame.


* ``laser.position`` (`3 floats in 3D and 2D` ; in meters)
    The coordinates of one of the point of the antenna that will emit the laser.
    The plane of the antenna is entirely defined by ``laser.position`` and ``laser.direction``.

    ``laser.position`` also corresponds to the origin of the coordinates system
    for the laser tranverse profile. For instance, for a Gaussian laser profile,
    the peak of intensity will be at the position given by ``laser.position``.
    This variable can thus be used to shift the position of the laser pulse
    transversally.

    .. note::
        In 2D, ``laser.position`` is still given by 3 numbers, but the second number is ignored.

    When running a **boosted-frame simulation**, provide the value of
    ``laser.position`` in the laboratory frame, and use ``warpx.gamma_boost``
    to automatically perform the conversion to the boosted frame. Note that,
    in this case, the laser antenna will be moving, in the boosted frame.

*  ``laser.profile_t_peak`` (`float`; in seconds)
    The time at which the laser reaches its peak intensity, at the position
    given by ``laser.position`` (only used for the ``"gaussian"`` profile)

    When running a **boosted-frame simulation**, provide the value of
    ``laser.profile_t_peak`` in the laboratory frame, and use ``warpx.gamma_boost``
    to automatically perform the conversion to the boosted frame.

*  ``laser.profile_duration`` (`float` ; in seconds)

    The duration of the laser, defined as :math:`\tau` below:

    - For the ``"gaussian"`` profile:

    .. math::

        E(\boldsymbol{x},t) \propto \exp\left( -\frac{(t-t_{peak})^2}{\tau^2} \right)

    - For the ``"harris"`` profile:

    .. math::

        E(\boldsymbol{x},t) \propto \frac{1}{32}\left[10 - 15 \cos\left(\frac{2\pi t}{\tau}\right) + 6 \cos\left(\frac{4\pi t}{\tau}\right) - \cos\left(\frac{6\pi t}{\tau}\right) \right]\Theta(\tau - t)

    When running a **boosted-frame simulation**, provide the value of
    ``laser.profile_duration`` in the laboratory frame, and use ``warpx.gamma_boost``
    to automatically perform the conversion to the boosted frame.

* ``laser.profile_waist`` (`float` ; in meters)
    The waist of the transverse Gaussian laser profile, defined as :math:`w_0` :

    .. math::

        E(\boldsymbol{x},t) \propto \exp\left( -\frac{\boldsymbol{x}_\perp^2}{w_0^2} \right)

* ``laser.wavelength`` (`float`; in meters)
    The wavelength of the laser in vacuum.

    When running a **boosted-frame simulation**, provide the value of
    ``laser.wavelength`` in the laboratory frame, and use ``warpx.gamma_boost``
    to automatically perform the conversion to the boosted frame.

* ``laser.polarization`` (`3 floats in 3D and 2D`)
    The coordinates of a vector that points in the direction of polarization of
    the laser. The norm of this vector is unimportant, only its direction matters.

    .. note::
        Even in 2D, all the 3 components of this vectors are important (i.e.
        the polarization can be orthogonal to the plane of the simulation).

*  ``laser.direction`` (`3 floats in 3D`)
    The coordinates of a vector that points in the propagation direction of
    the laser. The norm of this vector is unimportant, only its direction matters.

    The plane of the antenna that will emit the laser is orthogonal to this vector.

    .. warning::

        When running **boosted-frame simulations**, ``laser.direction`` should
        be parallel to ``warpx.boost_direction``, for now.

* ``laser.profile_focal_distance`` (`float`; in meters)
    The distance from ``laser_position`` to the focal plane.
    (where the distance is defined along the direction given by ``laser.direction``.)

    Use a negative number for a defocussing laser instead of a focussing laser.

    When running a **boosted-frame simulation**, provide the value of
    ``laser.profile_focal_distance`` in the laboratory frame, and use ``warpx.gamma_boost``
    to automatically perform the conversion to the boosted frame.


Numerics and algorithms
-----------------------

* ``warpx.use_filter`` (`0 or 1`)
    Whether to smooth the charge and currents on the mesh, after depositing
    them from the macroparticles. This uses a bilinear filter
    (see the sub-section **Filtering** in :doc:`../theory/theory`).

* ``algo.current_deposition`` (`integer`)
    The algorithm for current deposition:

     - ``0``: Esirkepov deposition, vectorized
     - ``1``: Esirkepov deposition, non-optimized
     - ``2``: Direct deposition, vectorized
     - ``3``: Direct deposition, non-optimized

     .. warning ::

        The vectorized Esirkepov deposition
        (``algo.current_deposition=0``) is currently not functional in WarpX.
        All the other methods (``1``, ``2`` and ``3``) are functional.

* ``algo.charge_deposition`` (`integer`)
    The algorithm for the charge density deposition:

     - ``0``: Vectorized version
     - ``1``: Non-optimized version

* ``algo.field_gathering`` (`integer`)
    The algorithm for field gathering:

     - ``0``: Vectorized version
     - ``1``: Non-optimized version

* ``algo.particle_pusher`` (`integer`)
    The algorithm for the particle pusher:

     - ``0``: Boris pusher
     - ``1``: Vay pusher

* ``interpolation.nox``, ``interpolation.noy``, ``interpolation.noz`` (`integer`)
    The order of the shape factors for the macroparticles, for the 3 dimensions of space. Lower-order shape factors result in faster simulations, but more noisy results,

    Note that the implementation in WarpX is more efficient when these 3 numbers are equal, and when they are between 1 and 3.

Diagnostics and output
----------------------

* ``amr.plot_int`` (`integer`)
    The number of PIC cycles inbetween two consecutive data dumps. Use a
    negative number to disable data dumping.

* ``warpx.do_boosted_frame_diagnostic`` (`0 or 1`)
    Whether to use the **back-transformed diagnostics** (i.e. diagnostics that
    perform on-the-fly conversion to the laboratory frame, when running
    boosted-frame simulations)

* ``warpx.num_snapshots_lab`` (`integer`)
    Only used when ``warpx.do_boosted_frame_diagnostic`` is ``1``.
    The number of lab-frame snapshots that will be written.

* ``warpx.dt_snapshots_lab`` (`float`, in seconds)
    Only used when ``warpx.do_boosted_frame_diagnostic`` is ``1``.
    The time interval inbetween the lab-frame snapshots (where this
    time interval is expressed in the laboratory frame).
