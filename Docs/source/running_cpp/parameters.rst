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

* ``warpx.boost_direction`` (string: ``x``, ``y`` or ``z``)
    The direction of the Lorentz-transform for boosted-frame simulations
    (The direction ``y`` cannot be used in 2D simulations.)

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

* ``warpx.fine_tag_lo`` and ``warpx.fine_tag_hi`` (`2 floats in 2D`, `3 integers in 3D`; in meters) optional
    **When using static mesh refinement with 1 level**, the extent of the refined patch.
    This patch is rectangular, and thus its extent is given here by the coordinates
    of the lower corner (``warpx.fine_tag_lo``) and upper corner (``warpx.fine_tag_hi``).

Distribution across MPI ranks and parallelization
-------------------------------------------------


* ``amr.max_grid_size`` (`integer`) optional (default `128`)
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

* ``warpx.load_balance_int`` (`integer`) optional (default `-1`)
    How often WarpX should try to redistribute the work across MPI ranks,
    in order to have better load balancing (expressed in number of PIC cycles
    inbetween two consecutive attempts at redistributing the work).
    Use 0 to disable load_balancing.

    When performing load balancing, WarpX measures the wall time for
    computational parts of the PIC cycle. It then uses this data to decide
    how to redistribute the subdomains across MPI ranks. (Each subdomain
    is unchanged, but its owner is changed in order to have better performance.)
    This relies on each MPI rank handling several (in fact many) subdomains
    (see ``max_grid_size``).

* ``warpx.load_balance_with_sfc`` (`0` or `1`) optional (default `0`)
    If this is `1`: use a Space-Filling Curve (SFC) algorithm in order to perform load-balancing of the simulation.
    If this is `0`: the Knapsack algorithm is used instead.

* ``warpx.do_dynamic_scheduling`` (`0` or `1`) optional (default `1`)
    Whether to activate OpenMP dynamic scheduling.

Math parser and user-defined constants
--------------------------------------

WarpX provides a math parser that reads expressions in the input file.
It can be used to define the plasma density profile, the plasma momentum
distribution or the laser field (see below `Particle initialization` and
`Laser initialization`).

The parser reads python-style expressions between double quotes, for instance
``"a0*x**2 * (1-y*1.e2) * (x>0)"`` is a valid expression where ``a0`` is a
user-defined constant and ``x`` and ``y`` are variables. The names are case sensitive. The factor
``(x>0)`` is `1` where `x>0` and `0` where `x<=0`. It allows the user to
define functions by intervals. User-defined constants can be used in parsed
functions only (i.e., ``density_function(x,y,z)`` and ``field_function(X,Y,t)``,
see below). User-defined constants can contain only letter, numbers and character _.
The name of each constant has to begin with a letter. The following names are used 
by WarpX, and cannot be used as user-defined constants: `x`, `y`, `z`, `X`, `Y`, `t`.
For example, parameters ``a0`` and ``z_plateau`` can be specified with:

* ``my_constants.a0 = 3.0``
* ``my_constants.z_plateau = 150.e-6``

Particle initialization
-----------------------

* ``particles.nspecies`` (`int`)
    The number of species that will be used in the simulation.

* ``particles.species_names`` (`strings`, separated by spaces)
    The name of each species. This is then used in the rest of the input deck ;
    in this documentation we use `<species_name>` as a placeholder.

* ``particles.use_fdtd_nci_corr`` (`0` or `1`) optional (default `0`)
    Whether to activate the FDTD Numerical Cherenkov Instability corrector.

* ``particles.rigid_injected_species`` (`strings`, separated by spaces)
    List of species injected using the rigid injection method. For species injected
    using this method, particles are translated along the `+z` axis with constant velocity
    as long as their ``z`` coordinate verifies ``z<zinject_plane``. When ``z>zinject_plane``,
    particles are pushed in a standard way, using the specified pusher.

* ``<species_name>.charge`` (`float`)
    The charge of one `physical` particle of this species.

* ``<species_name>.mass`` (`float`)
    The mass of one `physical` particle of this species.

* ``<species_name>.injection_style`` (`string`)
    Determines how the particles will be injected in the simulation.
    The options are:

    * ``NUniformPerCell``: injection with a fixed number of evenly-spaced particles per cell.
      This requires the additional parameter ``<species_name>.num_particles_per_cell_each_dim``.

    * ``NRandomPerCell``: injection with a fixed number of randomly-distributed particles per cell.
      This requires the additional parameter ``<species_name>.num_particles_per_cell``.

* ``<species_name>.profile`` (`string`)
    Density profile for this species. The options are:

    * ``constant``: Constant density profile within the box, or between ``<species_name>.xmin``
      and ``<species_name>.xmax`` (and same in all directions). This requires additional
      parameter ``<species_name>.density``. i.e., the plasma density in :math:`m^{-3}`.

    * ``parse_density_function``: the density is given by a function in the input file.
      It requires additional argument ``<species_name>.density_function(x,y,z)``, which is a
      mathematical expression for the density of the species, e.g.
      ``electrons.density_function(x,y,z) = "n0+n0*x**2*1.e12"`` where ``n0`` is a
      user-defined constant, see above.

* ``<species_name>.momentum_distribution_type`` (`string`)
    Distribution of the normalized momentum (`u=p/mc`) for this species. The options are:

    * ``constant``: constant momentum profile. This requires additional parameters
      ``<species_name>.ux``, ``<species_name>.uy`` and ``<species_name>.uz``, the normalized
      momenta in the x, y and z direction respectively.

    * ``gaussian``: gaussian momentum distribution in all 3 directions. This requires
      additional arguments for the average momenta along each direction
      ``<species_name>.ux_m``, ``<species_name>.uy_m`` and ``<species_name>.uz_m`` as
      well as standard deviations along each direction ``<species_name>.ux_th``,
      ``<species_name>.uy_th`` and ``<species_name>.uz_th``.

    * ``radial_expansion``: momentum depends on the radial coordinate linearly. This
      requires additional parameter ``u_over_r`` which is the slope.

    * ``parse_momentum_function``: the momentum is given by a function in the input
      file. It requires additional arguments ``<species_name>.momentum_function_ux(x,y,z)``,
      ``<species_name>.momentum_function_uy(x,y,z)`` and ``<species_name>.momentum_function_uz(x,y,z)``,
      which gives the distribution of each component of the momentum as a function of space.

* ``<species_name>.zinject_plane`` (`float`)
    Only read if  ``<species_name>`` is in ``particles.rigid_injected_species``.
    Injection plane when using the rigid injection method.
    See ``particles.rigid_injected_species`` above.

* ``<species_name>.rigid_avance`` (`bool`)
    Only read if ``<species_name>`` is in ``particles.rigid_injected_species``.

    * If ``false``, each particle is advanced with its
      own velocity ``vz`` until it reaches ``zinject_plane``.

    * If ``true``, each particle is advanced with the average speed of the species
      ``vzbar`` until it reaches ``zinject_plane``.

* ``<species_name>.do_backward_injection`` (`bool`)
    Inject a backward-propagating beam to reduce the effect of charge-separation
    fields when running in the boosted frame. See examples.

* ``<species_name>.do_splitting`` (`bool`) optional (default `0`)
    Split particles of the species when crossing the boundary from a lower
    resolution domain to a higher resolution domain.

* ``<species_name>.split_type`` (`int`) optional (default `0`)
    Splitting technique. When `0`, particles are split along the simulation
    axes (4 particles in 2D, 6 particles in 3D). When `1`, particles are split
    along the diagonals (4 particles in 2D, 8 particles in 3D).

* ``warpx.serialize_ics`` (`0 or 1`)
    Whether or not to use OpenMP threading for particle initialization.

Laser initialization
--------------------

* ``warpx.use_laser`` (`0 or 1`) optional (default `0`)
    Whether to activate the injection of a laser pulse in the simulation

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

* ``laser.e_max`` (`float` ; in V/m)
    Peak amplitude of the laser field.

    For a laser with a wavelength :math:`\lambda = 0.8\,\mu m`, the peak amplitude
    is related to :math:`a_0` by:

    .. math::

        E_{max} = a_0 \frac{2 \pi m_e c}{e\lambda} = a_0 \times (4.0 \cdot 10^{12} \;V.m^{-1})

    When running a **boosted-frame simulation**, provide the value of ``laser.e_max``
    in the laboratory frame, and use ``warpx.gamma_boost`` to automatically
    perform the conversion to the boosted frame.

* ``laser.wavelength`` (`float`; in meters)
    The wavelength of the laser in vacuum.

    When running a **boosted-frame simulation**, provide the value of
    ``laser.wavelength`` in the laboratory frame, and use ``warpx.gamma_boost``
    to automatically perform the conversion to the boosted frame.

* ``laser.profile`` (`string`)
    The spatio-temporal shape of the laser. The options that are currently
    implemented are:

    - ``"Gaussian"``: The transverse and longitudinal profiles are Gaussian.
    - ``"Harris"``: The transverse profile is Gaussian, but the longitudinal profile
      is given by the Harris function (see ``laser.profile_duration`` for more details)
    - ``"parse_field_function"``: the laser electric field is given by a function in the
      input file. It requires additional argument ``laser.field_function(X,Y,t)``, which
      is a mathematical expression , e.g.
      ``laser.field_function(X,Y,t) = "a0*X**2 * (X>0) * cos(omega0*t)"`` where
      ``a0`` and ``omega0`` are a user-defined constant, see above. The profile passed
      here is the full profile, not only the laser envelope. ``t`` is time and ``X``
      and ``Y`` are coordinates orthogonal to ``laser.direction`` (not necessarily the
      x and y coordinates of the simulation). All parameters above are required, but
      none of the parameters below are used when ``laser.parse_field_function=1``. Even
      though ``laser.wavelength`` and ``laser.e_max`` should be included in the laser
      function, they still have to be specified as they are used for numerical purposes.

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

* ``laser.profile_focal_distance`` (`float`; in meters)
    The distance from ``laser_position`` to the focal plane.
    (where the distance is defined along the direction given by ``laser.direction``.)

    Use a negative number for a defocussing laser instead of a focussing laser.

    When running a **boosted-frame simulation**, provide the value of
    ``laser.profile_focal_distance`` in the laboratory frame, and use ``warpx.gamma_boost``
    to automatically perform the conversion to the boosted frame.

* ``laser.stc_direction`` (`3 floats`) optional (default `1. 0. 0.`)
    Direction of laser spatio-temporal couplings.
  	See definition in Akturk et al., Opt Express, vol 12, no 19 (2014).

* ``laser.zeta`` (`float`; in meters.seconds) optional (default `0.`)
    Spatial chirp at focus in direction ``laser.stc_direction``. See definition in
    Akturk et al., Opt Express, vol 12, no 19 (2014).

* ``laser.beta`` (`float`; in seconds) optional (default `0.`)
    Angular dispersion (or angular chirp) at focus in direction ``laser.stc_direction``.
    See definition in Akturk et al., Opt Express, vol 12, no 19 (2014).

* ``laser.phi2`` (`float`; in seconds**2) optional (default `0.`)
    Temporal chirp at focus.
    See definition in Akturk et al., Opt Express, vol 12, no 19 (2014).

Numerics and algorithms
-----------------------

* ``warpx.cfl`` (`float`)
    The ratio between the actual timestep that is used in the simulation
    and the CFL limit. (e.g. for `warpx.cfl=1`, the timestep will be
    exactly equal to the CFL limit.)

* ``warpx.use_filter`` (`0 or 1`)
    Whether to smooth the charge and currents on the mesh, after depositing
    them from the macroparticles. This uses a bilinear filter
    (see the sub-section **Filtering** in :doc:`../theory/theory`).

* ``warpx.filter_npass_each_dir`` (`3 int`) optional (default `1 1 1`)
    Number of passes along each direction for the bilinear filter.
    In 2D simulations, only the first two values are read.

* ``algo.current_deposition`` (`integer`)
    The algorithm for current deposition:

     - ``0``: Esirkepov deposition, vectorized
     - ``1``: Esirkepov deposition, non-optimized
     - ``2``: Direct deposition, vectorized
     - ``3``: Direct deposition, non-optimized

    .. warning::

        On GPU, use ``algo.current_deposition=0`` for Esirkepov
	or ``3`` for direct deposition.

* ``algo.charge_deposition`` (`integer`)
    The algorithm for the charge density deposition:

     - ``0``: Vectorized version
     - ``1``: Non-optimized version

* ``algo.field_gathering`` (`integer`)
    The algorithm for field gathering:

     - ``0``: Vectorized version
     - ``1``: Non-optimized version

    .. warning::

        The vectorized version does not run on GPU. Use
		``algo.field_gather=1`` when running on GPU.

* ``algo.particle_pusher`` (`integer`)
    The algorithm for the particle pusher:

     - ``0``: Boris pusher
     - ``1``: Vay pusher

* ``algo.maxwell_fdtd_solver`` (`string`)
    The algorithm for the FDTD Maxwell field solver:

     - ``yee``: Yee FDTD solver
     - ``ckc``: Cole-Karkkainen solver with Cowan
       coefficients (see Cowan - PRST-AB 16, 041303 (2013))

* ``interpolation.nox``, ``interpolation.noy``, ``interpolation.noz`` (`integer`)
    The order of the shape factors for the macroparticles, for the 3 dimensions of space.
    Lower-order shape factors result in faster simulations, but more noisy results,

    Note that the implementation in WarpX is more efficient when these 3 numbers are equal,
    and when they are between 1 and 3.

* ``psatd.nox``, ``psatd.noy``, ``pstad.noz`` (`integer`) optional (default `16` for all)
    The order of accuracy of the spatial derivatives, when using the code compiled with a PSATD solver.

* ``psatd.ngroups_fft`` (`integer`)
    The number of MPI groups that are created for the FFT, when using the code compiled with a PSATD solver.
    The FFTs are global within one MPI group and use guard cell exchanges in between MPI groups.
    (If ``ngroups_fft`` is larger than the number of MPI ranks used,
    than the actual number of MPI ranks is used instead.)

* ``psatd.fftw_plan_measure`` (`0` or `1`)
    Defines whether the parameters of FFTW plans will be initialized by
    measuring and optimizing performance (``FFTW_MEASURE`` mode; activated by default here).
    If ``psatd.fftw_plan_measure`` is set to ``0``, then the best parameters of FFTW
    plans will simply be estimated (``FFTW_ESTIMATE`` mode).
    See `this section of the FFTW documentation <http://www.fftw.org/fftw3_doc/Planner-Flags.html>`__
    for more information.


Diagnostics and output
----------------------

* ``amr.plot_int`` (`integer`)
    The number of PIC cycles inbetween two consecutive data dumps. Use a
    negative number to disable data dumping.

* ``warpx.dump_plotfiles`` (`0` or `1`) optional
    Whether to dump the simulation data in
    `AMReX plotfile <https://amrex-codes.github.io/amrex/docs_html/IO.html>`__
    format. This is ``1`` by default, unless WarpX is compiled with openPMD support.

* ``warpx.dump_openpmd`` (`0` or `1`) optional
    Whether to dump the simulation data in
    `openPMD <https://github.com/openPMD>`__ format.
    When WarpX is compiled with openPMD support, this is ``1`` by default.

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

* ``warpx.plot_raw_fields`` (`0` or `1`) optional (default `0`)
    By default, the fields written in the plot files are averaged on the nodes.
    When ```warpx.plot_raw_fields`` is `1`, then the raw (i.e. unaveraged)
    fields are also saved in the plot files.

* ``warpx.plot_raw_fields_guards`` (`0` or `1`)
    Only used when ``warpx.plot_raw_fields`` is ``1``.
    Whether to include the guard cells in the output of the raw fields.

* ``warpx.plot_finepatch`` (`0` or `1`)
    Only used when mesh refinement is activated and ``warpx.plot_raw_fields`` is ``1``.
    Whether to output the data of the fine patch, in the plot files.

* ``warpx.plot_crsepatch`` (`0` or `1`)
    Only used when mesh refinement is activated and ``warpx.plot_raw_fields`` is ``1``.
    Whether to output the data of the coarse patch, in the plot files.

* ``warpx.plot_coarsening_ratio`` (`int` ; default: `1`)
    Reduce size of the field output by this ratio in each dimension.
    (This is done by averaging the field.) ``plot_coarsening_ratio`` should
    be an integer divisor of ``blocking_factor``.

* ``warpx.particle_plot_vars`` (`strings`, separated by spaces ; default: all)
    Control which particle variables get written to the plot file. Choices are:
    `w`, `ux`, `uy`, `uz`, `Ex`, `Ey`, `Ez`, `Bx`, `By`, and `Bz`.
    The particle positions and ids are always included.

* ``amr.plot_file`` (`string`)
    Root for output file names. Supports sub-directories. Default `diags/plotfiles/plt`

Checkpoints and restart
-----------------------
WarpX supports checkpoints/restart via AMReX.

* ``amr.check_int`` (`integer`)
    The number of iterations between two consecutive checkpoints. Use a
    negative number to disable checkpoints.

* ``amr.restart`` (`string`)
    Name of the checkpoint file to restart from. Returns an error if the folder does not exist
    or if it is not properly formatted.
