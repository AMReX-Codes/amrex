
.. _sec:basics:debugging:

Debugging
=========

Debugging is an art.  Everyone has their own favorite method.  Here we
offer a few tips we have found to be useful.

To help debugging, AMReX handles various signals in the C standard
library raised in the runs.  This gives us a chance to print out more
information using Linux/Unix backtrace capability.  The signals
include segmentation fault (or "segfault"), interruption by the user (control-c), assertion
errors, and floating point exceptions (NaNs, divided by zero and
overflow).  The handling of segfault, assertion errors and
interruption by control-C are enabled by default.  Note that
``AMREX_ASSERT()`` is only on when compiled with ``DEBUG=TRUE`` or
``USE_ASSERTION=TRUE`` in GNU make, or with ``-DCMAKE_BUILD_TYPE=Debug`` or
``-DAMReX_ASSERTIONS=YES`` in CMake.  The trapping of floating point exceptions is not
enabled by default unless the code is compiled with ``DEBUG=TRUE`` in GNU make, or with
``-DCMAKE_BUILD_TYPE=Debug`` or ``-DAMReX_FPE=YES`` in CMake to turn on compiler flags
if supported.  Alternatively, one can always use runtime parameters to control the
handling of floating point exceptions: ``amrex.fpe_trap_invalid`` for
NaNs, ``amrex.fpe_trap_zero`` for division by zero and
``amrex.fpe_trap_overflow`` for overflow.  To more effectively trap the
use of uninitialized values, AMReX also initializes ``FArrayBox``\ s in
``MulitFab``\ s and arrays allocated by ``bl_allocate`` to signaling NaNs when it is compiled
with ``TEST=TRUE`` or ``DEBUG=TRUE`` in GNU make, or with ``-DCMAKE_BUILD_TYPE=Debug`` in CMake.
One can also control the setting for ``FArrayBox`` using the runtime parameter, ``fab.init_snan``.

One can get more information than the backtrace of the call stack by
instrumenting the code.  Here is an example.
You know the line ``Real rho = state(cell,0);`` is causing a segfault.  You
could add a print statement before that.  But it might print out
thousands (or even millions) of line before it hits the segfault.  What
you could do is the following,

.. highlight:: c++

::

   #include <AMReX_BLBackTrace.H>

   std::ostringstream ss;
   ss << "state.box() = " << state.box() << " cell = " << cell;
   BL_BACKTRACE_PUSH(ss.str()); // PUSH takes std::string

   Real rho = state(cell,0);  // state is a Fab, and cell is an IntVect.

   BL_BACKTRACE_POP(); // One can omit this line.  In that case,
                       // there is an implicit POP when "PUSH" is
                       // out of scope.

When it hits the segfault, you will only see the last print out.

Writing a ``MultiFab`` to disk with

.. highlight:: c++

::

    VisMF::Write(const FabArray<FArrayBox>& mf, const std::string& name)

in ``AMReX_VisMF.H`` and examining it with ``Amrvis`` (section
:ref:`sec:amrvis`) can be helpful as well.  In
``AMReX_MultiFabUtil.H``, function

.. highlight:: c++

::

    void print_state(const MultiFab& mf, const IntVect& cell, const int n=-1,
                     const IntVect& ng = IntVect::TheZeroVector());

can output the data for a single cell. ``n`` is the component, with the default being
to print all components. ``ng`` is the number of ghost cells to include.

Valgrind is one of our favorite debugging tools.  For MPI runs, one can
tell Valgrind to output to different files for different processes.
For example,

.. highlight:: console

::

    mpiexec -n 4 valgrind --leak-check=yes --track-origins=yes --log-file=vallog.%p ./foo.exe ...

Breaking into Debuggers
-----------------------

In order to break into debuggers and use modern IDEs, the backtrace signal handling described above needs to be disabled.

The following runtime options need to be set in order to prevent AMReX from catching the break signals before a debugger can attach to a crashing process:

::

   amrex.throw_exception = 1
   amrex.signal_handling = 0

This default behavior can also be modified by applications, see for example `this custom application initializer <https://github.com/Exawind/amr-wind/blob/84f81a990152f4f748c1ab0fa17c8c663e51df86/amr-wind/main.cpp#L21>`__.


.. _sec:gpu:debugging:

Basic Gpu Debugging
-------------------


The asynchronous nature of GPU execution can make tracking down bugs complex.
The relative timing of improperly coded functions can cause variations in output and the timing of error messages
may not linearly relate to a place in the code.
One strategy to isolate specific kernel failures is to add ``amrex::Gpu::synchronize()`` or ``amrex::Gpu::streamSynchronize()`` after every ``ParallelFor`` or similar ``amrex::launch`` type call.
These synchronization commands will halt execution of the code until the GPU or GPU stream, respectively, has finished processing all previously requested tasks, thereby making it easier to locate and identify sources of error.

Debuggers and Related Tools
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Users may also find debuggers useful. Architecture agnostic tools include ``gdb``, ``hpctoolkit``, and ``Valgrind``. Note that there are architecture specific implementations of ``gdb`` such as ``cuda-gdb``, ``rocgdb``, ``gdb-amd``, and the Intel ``gdb``.
Usage of several of these variations are described in the following sections.

For advance debugging topics and tools, refer to system-specific documentation (e.g. https://docs.olcf.ornl.gov/systems/summit_user_guide.html#debugging).


CUDA-Specific Tests
^^^^^^^^^^^^^^^^^^^

- To test if your kernels have launched, run:

  ::

    nvprof ./main3d.xxx

  If using NVIDIA Nsight Compute instead, access ``nvprof`` functionality with:

  ::

    nsys nvprof ./main3d.xxx

- Run ``nvprof -o profile%p.nvvp ./main3d.xxxx`` or
  ``nsys profile -o nsys_out.%q{SLURM_PROCID}.%q{SLURM_JOBID} ./main3d.xxx`` for
  a small problem and examine page faults using ``nvvp`` or ``nsight-sys $(pwd)/nsys_out.#.######.qdrep``.

- Run under ``cuda-memcheck`` or the newer version ``compute-sanitizer`` to identify memory errors.

- Run under ``cuda-gdb`` to identify kernel errors.

- To help identify race conditions, globally disable asynchronicity of kernel launches for all
  CUDA applications by setting ``CUDA_LAUNCH_BLOCKING=1`` in your environment variables. This
  will ensure that only one CUDA kernel will run at a time.

AMD ROCm-Specific Tests
^^^^^^^^^^^^^^^^^^^^^^^

- To test if your kernels have launched, run:

  ::

    rocprof ./main3d.xxx

- Run ``rocprof  --hsa-trace --stats --timestamp on --roctx-trace ./main3d.xxxx`` for
  a small problem and examine tracing using ``chrome://tracing``.

- Run under ``rocgdb`` for source-level debugging.

- To help identify if there are race conditions, globally disable asynchronicity of kernel launches by setting ``CUDA_LAUNCH_BLOCKING=1`` or ``HIP_LAUNCH_BLOCKING=1``
  in your environment variables. This will ensure only one kernel will run at a time.
  See the `AMD ROCm docs' chicken bits section`_ for more debugging environment variables.

.. _`AMD ROCm docs' chicken bits section`: https://rocmdocs.amd.com/en/latest/Programming_Guides/HIP_Debugging.html#chicken-bits

Intel GPU Specific Tests
^^^^^^^^^^^^^^^^^^^^^^^^

- To test if your kernels have launched, run:

  ::

    ./ze_tracer ./main3d.xxx

- Run Intel Advisor,
  ``advisor --collect=survey ./main3d.xxx`` for
  a small problem with 1 MPI process and examine metrics.

- Run under ``gdb`` with the `Intel Distribution for GDB`_.

- To report back-end information, set ``ZE_DEBUG=1`` in your environment variables.

.. _`Intel Distribution for GDB`: https://software.intel.com/content/www/us/en/develop/tools/oneapi/components/distribution-for-gdb.html
