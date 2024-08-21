
.. _chap:inputs:

Runtime Parameters
==================

.. role:: cpp(code)
   :language: c++

This chapter contains a list of AMReX :cpp:`ParmParse` runtime parameters
and their **default** values. They can be set by either including them in an
inputs file, or specifying them at the command line, or passing a function
to :cpp:`amrex::Initialize` and the function adds parameters to AMReX's
:cpp:`ParmParse`'s parameter database. For more information on
:cpp:`ParmParse`, see :ref:`sec:basics:parmparse`.

.. important:: AMReX reserves the following prefixes in :cpp:`ParmParse`
               parameters: ``amr``, ``amrex``, ``blprofiler``, ``device``,
               ``DistributionMapping``, ``eb2``, ``fab``, ``fabarray``,
               ``geometry``, ``particles``, ``tiny_profiler``, and
               ``vismf``.

AMR
---

AMReX applications with AMR use either :cpp:`class AmrCore` or the more
specialized :cpp:`class Amr`. Since :cpp:`class Amr` is derived from
:cpp:`class AmrCore`, the parameters for the :cpp:`AmrCore` class also apply
to the :cpp:`Amr` class. Additionally, :cpp:`class AmrCore` is derived from
:cpp:`class AmrMesh`, so :cpp:`AmrMesh` member functions are also available
to :cpp:`AmrCore` and :cpp:`Amr`.

AmrCore Class
^^^^^^^^^^^^^

Below are a list of important :cpp:`ParmParse` parameters. However, AMReX
applications can choose to avoid them entirely by use this :cpp:`AMRCore`
constructor :cpp:`AmrCore(Geometry const& level_0_geom, AmrInfo const&
amr_info)`, where :cpp:`struct AmrInfo` contains all the information that
can be set via :cpp:`ParmParse`.

.. py:data:: amr.verbose
   :type: int
   :value: 0

   This controls the verbosity level of :cpp:`AmrCore` functions.

.. py:data:: amr.n_cell
   :type: int array
   :value: [none]

   This parameter is used only when ``n_cell`` is not provided as an
   argument to :cpp:`AmrCore` constructors. It specifies the number of cells
   in each dimension on Level 0.

.. py:data:: amr.max_level
   :type: int
   :value: [none]

   This parameter is used only when ``max_level`` is not provided as an
   argument to :cpp:`AmrCore` constructors. It specifies the maximum level
   of refinement allowed. Note that the total number of levels, including
   the base level 0, is ``max_level+1``.

.. py:data:: amr.ref_ratio
   :type: int array
   :value: 2 2 2 ... 2

   If the refinement ratio is not provided as an argument to :cpp:`AmrCore`
   constructors and :py:data:`amr.ref_ratio_vect` is not found in the
   :cpp:`ParmParse` database, this parameter will be used to set the
   refinement ratios between AMR levels. If there are more AMR levels than
   the size of the integer parameter array, the last integer will be used as
   the refinement ratio for the unspecified levels. For example, if
   ``max_level`` is 4 and the provided ``amr.ref_ratio`` parameter is ``2
   4``, the refinement ratios are 2, 4, 4 and 4, for levels 0/1, 1/2, 2/3
   and 3/4, respectively.

.. py:data:: amr.ref_ratio_vect
   :type: int array
   :value: [none]

   If the refinement ratio is not provided as an argument to :cpp:`AmrCore`
   constructors and :py:data:`amr.ref_ratio_vect` is found in the
   :cpp:`ParmParse` database, it will be used to set the refinement ratios
   between AMR levels. It's an error if the size of the integer array, if
   found, is less than ``max_level*AMREX_SPACEDIM``. The first
   ``AMREX_SPACEDIM`` numbers specify the refinement ratios in the
   ``AMREX_SPACEDIM`` dimensions between levels 0 and 1, the next
   ``AMREX_SPACEDIM`` numbers specify the ratios for levels 1 and 2, and so
   on.

.. py:data:: amr.max_grid_size
   :type: int array
   :value: [build dependent]

   This controls the maximum grid size on AMR levels, one value for each
   level. If the size of the integer array is less than the total number of
   levels, the last integer will be used for the unspecified levels. The
   default value is 128 for 1D and 2D runs. For 3D runs, the default value
   is 64 and 32, for GPU and CPU runs, respectively. Note that the user can
   also call :cpp:`AmrMesh::SetMaxGridSize` to set the maximum grid
   sizes. Additionally, the values set by this parameter can be overridden
   by :py:data:`amr.max_grid_size_x`, :py:data:`amr.max_grid_size_y` and
   :py:data:`amr.max_grid_size_z`.

.. py:data:: amr.max_grid_size_x
   :type: int array
   :value: [none]

   If provided, this will override the maximum grid size in the x-direction
   set by :py:data:`amr.max_grid_size`. If the size of the integer array is
   less than the total number of levels, the last integer will be used for
   the unspecified levels.

.. py:data:: amr.max_grid_size_y
   :type: int array
   :value: [none]

   If provided, this will override the maximum grid size in the y-direction
   set by :py:data:`amr.max_grid_size`. If the size of the integer array is
   less than the total number of levels, the last integer will be used for
   the unspecified levels.

.. py:data:: amr.max_grid_size_z
   :type: int array
   :value: [none]

   If provided, this will override the maximum grid size in the z-direction
   set by :py:data:`amr.max_grid_size`. If the size of the integer array is
   less than the total number of levels, the last integer will be used for
   the unspecified levels.

.. py:data:: amr.blocking_factor
   :type: int array
   :value: [build dependent]

   This controls the blocking factor on AMR levels, one value for each
   level. If the size of the integer array is less than the total number of
   levels, the last integer will be used for the unspecified levels. The
   default value is 8. Note that the user can also call
   :cpp:`AmrMesh::SetBlockingFactor` to set the blocking
   factors. Additionally, the values set by this parameter can be overridden
   by :py:data:`amr.blocking_factor_x`, :py:data:`amr.blocking_factor_y` and
   :py:data:`amr.blocking_factor_z`.

.. py:data:: amr.blocking_factor_x
   :type: int array
   :value: [none]

   If provided, this will override the blocking factor in the x-direction
   set by :py:data:`amr.blocking_factor`. If the size of the integer array
   is less than the total number of levels, the last integer will be used
   for the unspecified levels.

.. py:data:: amr.blocking_factor_y
   :type: int array
   :value: [none]

   If provided, this will override the blocking factor in the y-direction
   set by :py:data:`amr.blocking_factor`. If the size of the integer array
   is less than the total number of levels, the last integer will be used
   for the unspecified levels.

.. py:data:: amr.blocking_factor_z
   :type: int array
   :value: [none]

   If provided, this will override the blocking factor in the z-direction
   set by :py:data:`amr.blocking_factor`. If the size of the integer array
   is less than the total number of levels, the last integer will be used
   for the unspecified levels.

.. py:data:: amr.n_proper
   :type: int
   :value: 1

   This parameter controls the proper nesting of grids on AMR levels. For
   example, if we have ``blocking_factor = 8``, ``ref_ratio = 2`` and
   ``n_proper = 1``, there will be at least ``8/2*1 = 4`` coarse level cells
   outside the fine level grids except at the physical boundaries. Note that
   the user can also call :cpp:`AmrMesh::SetNProper(int)` to set the proper
   nesting parameter.

.. py:data:: amr.grid_eff
   :type: amrex::Real
   :value: 0.7

   This parameter controls the grid efficiency threshold during grid
   creation. While a higher value can enhance efficiency, it may negatively
   impact overall performance, especially for GPU runs, because it tends to
   create smaller grids. Note that the user can also call
   :cpp:`AmrMesh::SetGridEff(Real)` to set the grid efficiency threshold.

.. py:data:: amr.n_error_buf
   :type: int array
   :value: 1 1 1 ... 1

   This parameter controls how many extra cells will be tagged around every
   tagged cell. For example, if ``n_error_buf = 2``, tagging cell
   ``(i,j,k)`` will result in the tagging of the region of from lower corner
   ``(i-2,j-2,k-2)`` to upper corner ``(i+2,j+2,k+2)``. If the size of the
   integer array is less than the number of levels, the last integer will be
   used for the unspecified levels. Note that the values set by this
   parameter can be overridden by :py:data:`amr.n_error_buf_x`,
   :py:data:`amr.n_error_buf_y` and :py:data:`amr.n_error_buf_z`.


.. py:data:: amr.n_error_buf_x
   :type: int array
   :value: [none]

   This parameter controls the error buffer size in the x-direction. If the
   size of the integer array is less than the number of levels, the last
   integer will be used for the unspecified levels.

.. py:data:: amr.n_error_buf_y
   :type: int array
   :value: [none]

   This parameter controls the error buffer size in the y-direction. If the
   size of the integer array is less than the number of levels, the last
   integer will be used for the unspecified levels.

.. py:data:: amr.n_error_buf_z

   This parameter controls the error buffer size in the z-direction. If the
   size of the integer array is less than the number of levels, the last
   integer will be used for the unspecified levels.

.. py:data:: amr.refine_grid_layout
   :type: bool
   :value: true

   If it's true, AMReX will attempt to chop new grids into smaller chunks
   ensuring at least one grid per MPI process, provided this does not
   violate the blocking factor constraint.

.. py:data:: amr.refine_grid_layout_x
   :type: bool
   :value: [none]

   This parameter, if found, will override the
   :py:data:`amrex.refine_grid_layout` parameter in the x-direction.

.. py:data:: amr.refine_grid_layout_y
   :type: bool
   :value: [none]

   This parameter, if found, will override the
   :py:data:`amrex.refine_grid_layout` parameter in the y-direction.

.. py:data:: amr.refine_grid_layout_z
   :type: bool
   :value: [none]

   This parameter, if found, will override the
   :py:data:`amrex.refine_grid_layout` parameter in the z-direction.

.. py:data:: amr.check_input
   :type: bool
   :value: true

   If this is true, AMReX will check if the various parameters in
   :cpp:`AmrMesh` are reasonable.

Amr Class
^^^^^^^^^

.. warning:: These parameters are specific to :cpp:`class Amr` based
             applications. If your application use :cpp:`class AmrCore`
             directly, they do not apply unless you have provided
             implementations for them.

Subcycling
""""""""""

.. py:data:: amr.subcycling_mode
   :type: string
   :value: Auto

   This controls the subcycling mode of :cpp:`class Amr`. Possible value
   are ``None`` for no subcycling, or ``Auto`` for subcycling.

Regrid
""""""

.. py:data:: amr.regrid_int
   :type: int array
   :value: 1 1 1 ... 1

   This controls how often we perform the regrid operation on AMR levels 0
   to ``max_level-1``. If the parameter is a single value, it will be used
   on all levels. If the parameter is an array of more than one values, the
   size must be at least ``max_level`` and values after the first
   ``max_level`` elements are ignored.

.. py:data:: amr.regrid_on_restart
   :type: bool
   :value: false

   This controls whether we perform regrid immediately after restart.

.. py:data:: amr.force_regrid_level_zero
   :type: bool
   :value: false

   This controls whether we perform regrid on level 0.

.. py:data:: amr.compute_new_dt_on_regrid
   :type: bool
   :value: false

   This controls whether we re-compute ``dt`` after regrid.

.. py:data:: amr.initial_grid_file
   :type: string
   :value: [none]

   If this is set, the initial grids will be read from the specified file.

.. py:data:: amr.regrid_file
   :type: string
   :value: [none]

   If this is set, regrid will use the grids in the specified file.

I/O
"""

.. py:data:: amr.restart
   :type: string
   :value: [none]

   If this is set, the simulation will restart from the specified checkpoint
   file.

.. py:data:: amr.plotfile_on_restart
   :type: bool
   :value: false

   If this is set to true, a plotfile will be written after restart.

.. py:data:: amr.file_name_digits
   :type: int
   :value: 5

   This parameter specifies the minimum number of digits in checkpoint and
   plotfile names.

.. py:data:: amr.checkpoint_files_output
   :type: bool
   :value: true

   This controls whether we write checkpoint files.

.. py:data:: amr.check_file
   :type: string
   :value: chk

   This sets the "root" of checkpoint file names. For example, the
   checkpoint files are named ``chk00000``, ``chk001000``, etc. by default.

.. py:data:: amr.check_int
   :type: int
   :value: -1

   This controls the interval of writing checkpoint files, defined as the
   number of level 0 steps between each checkpoint. A value less than 1
   indicates no checkpoint files will be written.

.. py:data:: amr.check_per
   :type: amrex::Real
   :value: -1

   This controls the interval of writing checkpoint files, defined as the
   time (not the wall time) elapsed between each checkpoint. A value less
   or equal to 0 indicates no checkpoint files will be written.

.. py:data:: amr.checkpoint_nfiles
   :type: int
   :value: 64

   This is the maximum number of binary files per :cpp:`MultiFab` when
   writing checkpoint files.

.. py:data:: amr.plot_files_output
   :type: bool
   :value: true

   This controls whether we write plot files.

.. py:data:: amr.plot_file
   :type: string
   :value: plt

   This sets the "root" of plot file names. For example, the plot files are
   named ``plt00000``, ``plt001000``, etc. by default.

.. py:data:: amr.plot_int
   :type: int
   :value: -1

   This controls the interval of writing plot files, defined as the number
   of level 0 steps between each plot file. A value less than 1 indicates no
   plot files will be written.

.. py:data:: amr.plot_per
   :type: amrex::Real
   :value: -1

   This controls the interval of writing plot files, defined as the time
   (not the wall time) elapsed between each plot file. A value less or equal
   to 0 indicates no plot files will be written.

.. py:data:: amr.plot_log_per
   :type: amrex::Real
   :value: -1

   This controls the interval of writing plot files, defined as the
   ``log10`` time (not the wall time) elapsed between each plot file. A
   value less or equal to 0 indicates no plot files will be written.

.. py:data:: amr.plot_max_level
   :type: int
   :value: amr.max_level

   This controls the finest level in a plot file. For example, if the finest
   level in a run is 3, but this parameter is set to 1, only levels 0 and 1
   will be saved in a plot file.

.. py:data:: amr.plot_nfiles
   :type: int
   :value: 64

   This is the maximum number of binary files per :cpp:`MultiFab` when
   writing plot files.

.. py:data:: amr.plot_vars
   :type: string array
   :value: [none]

   If this parameter is set, the variables specified in the string array
   will be the state variables saved in the plot files. The special values
   ``ALL`` and ``NONE`` mean that all or none of the state variables will be
   saved. If this parameter is not set, all state variables will be saved.

.. py:data:: amr.derive_plot_vars
   :type: string array
   :value: [none]

   If this parameter is set, the variables specified in the string array
   will be the derive variables saved in the plot files. The special values
   ``ALL`` and ``NONE`` mean that all or none of the derive variables will
   be saved. If this parameter is not set, none of the derive variables will
   be saved.

.. py:data:: amr.small_plot_file
   :type: string
   :value: smallplt

   This sets the "root" of small plot file names. For example, the small
   plot files are named ``smallplt00000``, ``smallplt001000``, etc. by
   default.

.. py:data:: amr.small_plot_int
   :type: int
   :value: -1

   This controls the interval of writing small plot files, defined as the
   number of level 0 steps between each small plot file. A value less than 1
   indicates no small plot files will be written.

.. py:data:: amr.small_plot_per
   :type: amrex::Real
   :value: -1

   This controls the interval of writing small plot files, defined as the
   time (not the wall time) elapsed between each small plot file. A value
   less or equal to 0 indicates no small plot files will be written.

.. py:data:: amr.small_plot_log_per
   :type: amrex::Real
   :value: -1

   This controls the interval of writing small plot files, defined as the
   ``log10`` time (not the wall time) elapsed between each small plot
   file. A value less or equal to 0 indicates no small plot files will be
   written.

.. py:data:: amr.small_plot_vars
   :type: string array
   :value: [none]

   If this parameter is set, the variables specified in the string array
   will be the state variables saved in the small plot files. The special
   values ``ALL`` and ``NONE`` mean that all or none of the state variables
   will be saved. If this parameter is not set, none of the state variables
   will be saved.

.. py:data:: amr.derive_small_plot_vars
   :type: string array
   :value: [none]

   If this parameter is set, the variables specified in the string array
   will be the derive variables saved in the small plot files. The special
   values ``ALL`` and ``NONE`` mean that all or none of the derive variables
   will be saved. If this parameter is not set, none of the derive variables
   will be saved.

.. py:data:: amr.message_int
   :type: int
   :value: 10

   This controls the interval of checking messages during a run, defined as
   the number of level 0 steps between checks. A value less than 1 indicates
   no checking will be performed. A message refers to a file created by the
   user on the disk, where only the file name is checked, not its
   content. If the file name matches one of the following predefined names,
   appropriate actions will be taken.

   dump_and_continue
      Make a checkpoint file and continue running the simulation.

   stop_run
      Stop the simulation.

   dump_and_stop
      Make a checkpoint file and stop the simulation.

   plot_and_continue
      Make a plot file and continue running the simulation.

   small_plot_and_continue
      Make a small plot file and continue running the simulation.

.. py:data:: amr.write_plotfile_with_checkpoint
   :type: bool
   :value: true

   This parameter is for the message action discussed in
   :py:data:`amr.message_int`. It controls whether an action will make a
   plot file as well when asked to make a checkpoint file.

.. py:data:: amr.run_log
   :type: string
   :value: [none]

   If this parameter is set, the run log will be enabled and this is the log
   file name.

.. py:data:: amr.run_log_terse
   :type: string
   :value: [none]

   If this parameter is set, the terse run log will be enabled and this is
   the log file name.

.. py:data:: amr.grid_log
   :type: string
   :value: [none]

   If this parameter is set, the grid log will be enabled and this is the
   log file name.

.. py:data:: amr.data_log
   :type: string
   :value: [none]

   If this parameter is set, the data log will be enabled and this is the
   log file name.

Basic Controls
--------------

.. py:data:: amrex.verbose
   :type: int
   :value: 1

   This controls the verbosity level of AMReX. Besides using
   :cpp:`ParmParse`, you can also call :cpp:`amrex::SetVerbose(int)` to set
   it.

.. py:data:: amrex.init_snan
   :type: bool
   :value: [build dependent]

   This controls whether :cpp:`MultiFab`, :cpp:`FArrayBox`,
   :cpp:`BaseFab<double|float>`, :cpp:`PODVectors<double|float>`,
   :cpp:`Gpu::DeviceVector<double|float>`, etc. will be initialized to
   signaling NaNs at construction. The default value is true for debug
   builds. For non-debug builds, the default is false unless ``TEST=TRUE``
   for GNU Make or ``AMReX_TESTING`` is enabled for CMake.

.. py:data:: amrex.abort_on_unused_inputs
   :type: bool
   :value: false

   If this is true and there are unused :cpp:`ParmParse` parameters, AMReX
   will abort during :cpp:`amrex::Finalize`.

.. py:data:: amrex.parmparse.verbose
   :type: int
   :value: amrex.verbose

   If this is greater than zero, unused :cpp:`ParmParse` variables will be
   printed out during :cpp:`amrex::Finalize` or
   :cpp:`ParmParse::QueryUnusedInputs`. The parameter can also be set by
   calling :cpp:`amrex::ParmParse::SetVerbose(int)`.

.. py:data:: amrex.device.verbose
   :type: int
   :value: 0

   This controls whether AMReX prints out GPU device properties such name,
   vendor, total memory size, etc. This is only relevant for GPU runs.

.. py:data:: amrex.max_gpu_streams
   :type: int
   :value: 4

   This controls the number of GPU streams used by AMReX. It's only relevant
   for GPU runs.

.. py:data:: amrex.omp_threads
   :type: string
   :value: system

   If OpenMP is enabled, this can be used to set the default number of
   threads. Possible values are ``system``, ``nosmt``, or an integer
   string. The special value ``nosmt`` can be used to avoid using threads
   for virtual cores (aka Hyperthreading or SMT), as is default in OpenMP,
   and instead only spawns threads equal to the number of physical cores in
   the system.  For the values ``system`` and ``nosmt``, the environment
   variable ``OMP_NUM_THREADS`` takes precedence. If the string can be
   converted to an integer, ``OMP_NUM_THREADS`` is ignored.

.. py:data:: amrex.memory_log
   :type: string
   :value: memlog

   This is the name of the memory log file when memory profiling is enabled.

Communication
-------------

.. py:data:: amrex.use_gpu_aware_mpi
   :type: bool
   :value: false

   For GPU runs, this controls the memory type used for AMReX's
   communication buffers. When this is true, AMReX uses GPU device memory
   for communication data in MPI function calls. When this is false, the
   data are placed in pinned memory. Note that this flag does not enable
   GPU-aware MPI by itself. Enabling GPU-aware MPI is system
   dependent. Users should consult their system's documentation for
   instructions on setting up the environment and linking to GPU-aware MPI
   libraries.

Distribution Mapping
--------------------

.. py:data:: DistributionMapping.verbose
   :type: int
   :value: 0

   This controls the verbosity level of :cpp:`DistributionMapping`
   functions.

.. py:data:: DistributionMapping.strategy
   :type: string
   :value: SFC

   This is the default :cpp:`DistributionMapping` strategy. Possible values
   are ``SFC``, ``KNAPSACK``, ``ROUNDROBIN``, or ``RRSFC``. Note that the
   default strategy can also be set by calling
   :cpp:`DistributionMapping::strategy(DistributionMapping::Strategy)`.

Embedded Boundary
-----------------

.. py:data:: eb2.max_grid_size
   :type: int
   :value: 64

   This parameter specifies the maximum grid size in AMReX's internal EB
   database, not the user's data.

.. py:data:: eb2.extend_domain_face
   :type: bool
   :value: true

   This controls the behavior of the embedded boundary outside the
   domain. If this is true, the embedded boundary outside the domain is
   extended perpendicularly from the domain face. Otherwise, it's generated
   with the user provided implicit function. Note that this parameter can be
   overridden by the user when calling :cpp:`amrex::EB2::Build` with the
   optional parameter ``bool extend_domain_face``.

.. py:data:: eb2.num_coarsen_opt
   :type: int
   :value: 0

   If it is greater than 0, this parameter can speed up the EB
   generation. It indicates that the search for EB can be performed on grids
   coarsened by this factor and then the EB information details will be
   generated on the original grids. However, the user should be aware that
   setting this parameter too high could result in erroneous results. Also
   note that this parameter can be overridden by the user when calling
   :cpp:`amrex::EB2::Build` with the optional parameter ``int
   num_coarsen_opt``.

.. py:data:: eb2.geom_type
   :type: string
   :value: [none]

   There are two versions of the `amrex::EB2::Build` function that can be
   used to build EB. One version is a function template that takes a user
   provided :cpp:`GeometryShop`, while the other uses :cpp:`ParmParse`
   parameters to build EB. For the latter version, this parameter specifies
   the type of the EB. Possible values include the following.

   all_regular
      The entire domain is regular without any EB objects.

   parser
      The embedded boundary is describe by :py:data:`eb2.parser_function`.

   stl
      The embedded boundary will be built using an STL file specified by
      :py:data:`eb2.stl_file`.

.. py:data:: eb2.parser_function
   :type: string
   :value: [none]

   When ``eb2.geom_type = parser``, this parameter is a parser function
   string that contains a math expression describing the surface of the EB.

   .. seealso:: Section :ref:`sec:basics:parser`.

.. py:data:: eb2.stl_file
   :type: string
   :value: [none]

   When ``eb2.geom_type = stl``, this is a required string parameter
   specifying the STL file name.

.. py:data:: eb2.stl_scale
   :type: amrex:Real
   :value: 1

   When building EB using STL, the triangles in the STL file will be scaled
   by the given value of this optional parameter.

.. py:data:: eb2.stl_center
   :type: amrex::Real array
   :value: 0 0 0

   When building EB using STL, this optional parameter specifies the shifted
   center. The original coordinates in the STL file will be shifted by the
   provided values.

.. py:data:: eb2.stl_reverse_normal
   :type: bool
   :value: false

   When building EB using STL, the normal direction of the triangles in the
   STL file will be reversed if this optional parameter is set to true.

.. py:data:: eb2.small_volfrac
   :type: amrex::Real
   :value: [depend on the type of amrex::Real]

   This parameter specifies the threshold for small cells that will be
   converted to covered cells. The default value is ``1.e-14`` if
   :cpp:`amrex::Real` is ``double``, or ``1.e-5`` if :cpp:`amrex::Real` is
   ``float``.

.. py:data:: eb2.cover_multiple_cuts
   :type: bool
   :value: false

   If this parameter is set to true, multi-cut cells will be converted to
   covered cells.

   .. tip::  Because AMReX currently does not support multi-cut cells, it
             would be a runtime error if multi-cut cells are left unfixed.

.. py:data:: eb2.maxiter
   :type: int
   :value: 32

   Fixing small and multi-cut cells is an iterative process. This parameter
   specifies the maximum number of iterations for the fix-up process.

Error Handling
--------------

By default AMReX installs a signal handler that will be run when a signal
such as segfault is received. You can also enable floating point exception
trapping. The signal handler will print out backtraces that can be useful
for debugging.

.. note:: Floating point exception trapping is not enabled by default,
   because compilers might generate optimized SIMD code that raises the
   exceptions.

.. py:data:: amrex.signal_handling
   :type: bool
   :value: true

   This controls whether AMReX should handle signals.

.. py:data:: amrex.handle_sigsegv
   :type: bool
   :value: true

   If both this flag and ``amrex.signal_handling`` are true, ``SIGSEGV``
   will be handled by AMReX.

.. py:data:: amrex.handle_sigterm
   :type: bool
   :value: false

   If both this flag and ``amrex.signal_handling`` are true, ``SIGTERM``
   will be handled by AMReX. This flag is false by default because this
   could generate lots of backtrace files on some batch systems that issue
   ``SIGTERM`` for jobs running out of wall clock time.

.. py:data:: amrex.handle_sigint
   :type: bool
   :value: true

   If both this flag and ``amrex.signal_handling`` are true, ``SIGINT``
   will be handled by AMReX.

.. py:data:: amrex.handle_sigabrt
   :type: bool
   :value: true

   If both this flag and ``amrex.signal_handling`` are true, ``SIGABGT``
   will be handled by AMReX.

.. py:data:: amrex.handle_sigfpe
   :type: bool
   :value: true

   If both this flag and ``amrex.signal_handling`` are true, ``SIGFPE``
   will be handled by AMReX.

   .. seealso::
      Use :py:data:`amrex.fpe_trap_invalid`, :py:data:`amrex.fpe_trap_zero`
      and :py:data:`amrex.fpe_trap_overflow` to enable ``FE_INVALID``,
      ``FE_DIVBYZERO`` and ``FE_OVERFLOW`` trapping, respectively.

.. py:data:: amrex.handle_sigill
   :type: bool
   :value: true

   If both this flag and ``amrex.signal_handling`` are true, ``SIGILL``
   will be handled by AMReX.

.. py:data:: amrex.throw_exception
   :type: bool
   :value: false

   If this flag is true and ``amrex.signal_handling`` is false,
   :cpp:`amrex::Abort` and :cpp:`amrex::Error` will throw
   :cpp:`std::runtime_error` instead of aborting immediately. Note that
   according the C++ standard, if an exception is thrown and not caught,
   :cpp:`std::terminate` will be called.

.. py:data:: amrex.fpe_trap_invalid
   :type: bool
   :value: false

    If ``SIGFPE`` is handled by AMReX and this flag is true, ``FE_INVALID``
    (e.g., ``0/0``) trapping will be enabled. This flag has no effect on
    Windows.

.. py:data:: amrex.fpe_trap_zero
   :type: bool
   :value: false

    If ``SIGFPE`` is handled by AMReX and this flag is true,
    ``FE_DIVBYZERO`` (e.g., ``1/0``) trapping will be enabled. This flag has
    no effect on Windows.

.. py:data:: amrex.fpe_trap_overflow
   :type: bool
   :value: false

    If ``SIGFPE`` is handled by AMReX and this flag is true, ``FE_OVERFLOW``
    (i.e., the result is too large to be representable) trapping will be
    enabled. This flag has no effect on Windows.

Extern
------

Hypre
^^^^^

These parameters are relevant only when Hypre support is enabled.

.. py:data:: amrex.init_hypre
   :type: bool
   :value: true

   This controls whether AMReX should call ``HYPRE_Init()`` during
   :cpp:`amrex::Initialize`.

.. py:data:: amrex.hypre_spgemm_use_vendor
   :type: bool
   :value: false

   This controls whether HYPRE should use the vendor's ``SpGemm``
   functionality.

.. py:data:: amrex.hypre_spmv_use_vendor
   :type: bool
   :value: false

   This controls whether HYPRE should use the vendor's ``SpMV``
   functionality.

.. py:data:: amrex.hypre_sptrans_use_vendor
   :type: bool
   :value: false

   This controls whether HYPRE should use the vendor's ``SpTrans``
   functionality.

.. _sec:inputs:geom:

Geometry
--------

All these parameters are optional for constructing a :ref:`Geometry <sec:basics:geom>`
object. There are only used if the information is not provided via function
arguments.

.. py:data:: geometry.coord_sys
   :type: int
   :value: 0

   This specifies the coordinate system type with valid values being 0
   (Cartesian), or 1 (cylindrical), or 2 (spherical).

.. py:data:: geometry.prob_lo
   :type: amrex::Real array
   :value: 0 0 0

   This specifies the position of the lower corner of the physical domain.

.. py:data:: geometry.prob_hi
   :type: amrex::Real array
   :value: [none]

   This specifies the position of the upper corner of the physical
   domain. If this is provided, :py:data:`geometry.prob_extent` will be
   ignored.

.. py:data:: geometry.prob_extent
   :type: amrex::Real array
   :value: [none]

   This specifies the length of the physical domain. If
   :py:data:`geometry.prob_hi` is provided, this will be ignored.

.. py:data:: geometry.is_periodic
   :type: int array
   :value: 0 0 0

   These integer parameters are boolean flags to indicate whether the domain
   is periodic in each direction. It's considered true (i.e., periodic) if
   its value is non-zero, and false (i.e., non-periodic) if its value is
   zero.

I/O
---

.. py:data:: amrex.async_out
   :type: bool
   :value: false

   If this is true, AMReX's native mesh and particle plotfiles will be
   written asynchronously by a background thread.

.. py:data:: amrex.async_out_nfiles
   :type: into
   :value: 64

   This is the maximum number of binary files on each AMR level that will be
   used when AMReX writes a plotfile asynchronously.

.. py:data:: vismf.verbose
   :type: int
   :value: 0

   This controls the verbosity level of :cpp:`VisMF` functions.

Memory
------

.. py:data:: amrex.the_arena_init_size
   :type: long
   :value: [system dependent]

   This controls the main memory arena's initial size in bytes. For CPU
   runs, the default is 0, whereas for GPU runs, the default is set at run
   time to 3/4 of the system's device memory.

   .. tip:: Since ``amrex v24.08``, instead of
            ``amrex.the_arena_init_size=10000000000``, one can use
            ``amrex.the_arena_init_size=10'000'000'000`` or
            ``amrex.the_arena_init_size=1e10`` to set :cpp:`ParmParse`
            integer parameters like this one.

.. py:data:: amrex.the_device_arena_init_size
   :type: long
   :value: 8388608 [8 MB]

   This controls the GPU device arena's initial size in bytes. For CPU runs,
   this is ignored. If the main arena uses the device memory (as opposed to
   managed memory), this parameter is also ignored.

.. py:data:: amrex.the_managed_arena_init_size
   :type: long
   :value: 8388608 [8 MB]

   This controls the managed device arena's initial size in bytes. For CPU
   runs, this is ignored. If the main arena uses the managed memory (as
   opposed to device memory), this parameter is also ignored.

.. py:data:: amrex.the_pinned_arena_init_size
   :type: long
   :value: [system dependent]

   This controls the pinned host memory arena's initial size in bytes. The
   default is 8 MB for CPU runs. For GPU runs it's set to half of the GPU
   device memory by default.

.. py:data:: amrex.the_comms_arena_init_size
   :type: long
   :value: 8388608 [8 MB]

   This controls the MPI communication memory arena's initial size in bytes.

.. py:data:: amrex.the_arena_release_threshold
   :type: long
   :value: LONG_MAX

   This controls the release threshold of the main arena.

.. py:data:: amrex.the_device_arena_release_threshold
   :type: long
   :value: LONG_MAX

   This controls the release threshold of the device arena.

.. py:data:: amrex.the_managed_arena_release_threshold
   :type: long
   :value: LONG_MAX

   This controls the release threshold of the managed arena.

.. py:data:: amrex.the_pinned_arena_release_threshold
   :type: long
   :value: LONG_MAX

   This controls the release threshold of the pinned arena.

.. py:data:: amrex.the_comms_arena_release_threshold
   :type: long
   :value: LONG_MAX

   This controls the release threshold of the communication arena.

.. py:data:: amrex.the_async_arena_release_threshold
   :type: long
   :value: LONG_MAX

   This controls the release threshold of the asynchronous arena. Note that
   this is only relevant for the CUDA (>= 11.2) and HIP backends that
   support stream-ordered memory allocator.

.. py:data:: amrex.the_arena_is_managed
   :type: bool
   :value: false

   This controls if AMReX uses the managed memory for the main arena. This
   is only relevant for GPU runs.

.. py:data:: amrex.abort_on_out_of_gpu_memory
   :type: bool
   :value: false

   This controls if AMReX should simply abort when the reported free device
   memory is less than the amount an arena is asked to allocate. Note that
   for managed memory it's possible to allocate more than the amount of free
   device memory available. However, the code will be very slow. This
   parameter is only relevant for GPU runs.

.. py:data:: amrex.mf.alloc_single_chunk
   :type: bool
   :value: false

   This controls if all the data in a :cpp:`FabArray` (including
   :cpp:`MultiFab`) are in a contiguous chunk of memory.

.. py:data:: amrex.vector_growth_factor
   :type: amrex::Real
   :value: 1.5

   This controls the growth factor of :cpp:`amrex::PODVector` and its
   derived classes such as :cpp:`amrex::Gpu::DeviceVector`,
   :cpp:`amrex::Gpu::ManagedVector`, etc. A smaller value can avoid wasting
   memory, but it may result in a performance penalty during resizing.

Particles
---------

.. py:data:: particles.do_tiling
   :type: bool
   :value: false

   This controls whether tiling is enabled for particle containers.

.. py:data:: particles.tile_size
   :type: int array
   :value: 1024000 8 8

   When tiling is enabled, this is the default tile size. Note that a big
   number like 1024000 effectively turns tiling off in that direction.

.. py:data:: particles.do_mem_efficient_sort
   :type: bool
   :value: true

   This parameter controls whether the more memory efficient method will be
   used for sorting particles.

.. py:data:: particles.particles_nfiles
   :type: int
   :value: 256

   This is the maximum number of binary files per level for a particle
   container when writing checkpoint and plot files for particles. The
   special value of ``-1`` indicates one file per process.

Tiling
------

.. py:data:: fabarray.mfiter_tile_size
   :type: int array
   :value: [build dependent]

   This is the default size for :ref:`tiling <sec:basics:mfiter>`. For GPU
   runs, it is disabled by default. For CPU runs, it is disabled by default
   in 1D and 2D, but enabled in 3D with a tile size of 8 in the y and
   z-directions.

.. py:data:: fabarray.comm_tile_size
   :type: int array
   :value: [build dependent]

   This is the default tiling size used in moving data in and out of the MPI
   communication buffer . It is disabled by default for GPU runs, but
   enabled for CPU runs with a tile size of 8 in the y and z-directions (if
   they exist).

Tiny Profiler
-------------

These parameters are ignored unless profiling with :cpp:`TinyProfiler` is
enabled.

.. py:data:: tiny_profiler.verbose
   :type: int
   :value: 0

   If this value is greater than 0, messages about entering or leaving
   profiled regions will be printed on the I/O process.

.. py:data:: tiny_profiler.print_threshold
  :type: double
  :value: 1.0

  In the profiling report, regions with very small run times are not listed
  individually. Instead, they are included in a section named "Other". This
  parameter specifies the maximum inclusive run time that the "Other"
  section can take in percent relative to the total run time.

.. py:data:: tiny_profiler.device_synchronize_around_region
  :type: bool
  :value: false

  This parameter is only relevant for GPU runs. If it is set to true, the
  current GPU stream is synchronized when entering and leaving a profiling
  region. Because GPU kernels are asynchronous, time measurements without
  synchronization could be misleading. Enabling this parameter can provide
  more accurate measurements. However, the added synchronization points,
  which are unnecessary for correctness, could potentially degrade the
  performance.

.. py:data:: tiny_profiler.enabled
   :type: bool
   :value: true

   .. versionadded:: 24.09
      Runtime parameter `tiny_profiler.enabled``.

   This parameter can be used to disable tiny profiling including
   :cpp:`CArena` memory profiling at run time.

.. py:data:: tiny_profiler.memprof_enabled
   :type: bool
   :value: true

   .. versionadded:: 24.09
      Runtime parameter ``tiny_profiler.memprof_enabled``.

   This parameter can be used to disable :cpp:`CArena` memory profiling at
   run time. If ``tiny_profiler.enabled`` is false, this parameter has no
   effects.

.. py:data:: tiny_profiler.output_file
   :type: string
   :value: [empty]

   .. versionadded:: 24.09
      Runtime parameter ``tiny_profiler.output_file``.

   If this parameter is empty, the output of tiny profiling is dumped on the
   default out stream of AMReX. If it's not empty, it specifies the file
   name for the output. Note that ``/dev/null`` is a special name that means
   no output.
