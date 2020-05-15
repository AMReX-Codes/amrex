.. role:: cpp(code)
   :language: c++

.. role:: fortran(code)
   :language: fortran

Types of Profiling
==================

Currently you have two options for AMReX-specific profiling:
:ref:`sec:tiny:profiling` and :ref:`sec:full:profiling`.

.. _sec:tiny:profiling:

Tiny Profiling
----------------------

To enable "Tiny Profiling", if using GNU Make then set

::

  TINY_PROFILE = TRUE
  PROFILE      = FALSE

in your GNUMakefile.   If using cmake then set the following cmake flags

::

  AMREX_ENABLE_TINY_PROFILE = ON
  AMREX_ENABLE_BASE_PROFILE = OFF

Note that if you set ``PROFILE = TRUE``  (or ``AMREX_ENABLE_BASE_PROFILE =
ON``) then this will override the ``TINY_PROFILE`` flag and tiny profiling will
be disabled.

At the end of the run, a summary of exclusive and inclusive function times will
be written to stdout.  This output includes the minimum and maximum (over
processes) time spent in each routine as well as the average and the maximum
percentage of total run time.   See :ref:`sec:sample:tiny` for sample output.

The tiny profiler automatically writes the results to stdout at the end of your
code, when ``amrex::Finalize();`` is reached. However, you may want to write
partial profiling results to ensure your information is saved when you may fail
to converge or if you expect to run out of allocated time. Partial results can
be written at user-defined times by inserting the line:

::

  BL_PROFILE_TINY_FLUSH();

Any timers that have not reached their ``BL_PROFILE_VAR_STOP`` call or exited
their scope and deconstructed will not be included in these partial outputs.
(e.g., a properly instrumented ``main()`` should show a time of zero in all
partial outputs.) Therefore, it is recommended to place these flush calls in
easily identifiable regions of your code and outside of as many profiling
timers as possible, such as immediately before or after writing a checkpoint.

Also, since flush calls will print multiple, similar looking outputs to stdout,
it is also recommended to wrap any ``BL_PROFILE_TINY_FLUSH();`` calls in
informative ``amrex::Print()`` lines to ensure accurate identification of each
set of timers.

.. _sec:full:profiling:

Full Profiling
--------------

If you set ``PROFILE = TRUE`` then a ``bl_prof`` directory will be written that
contains detailed per-task timings for each processor.  This will be written in
``nfiles`` files (where ``nfiles`` is specified by the user).  In addition, an
exclusive-only set of function timings will be written to stdout.

Trace Profiling
~~~~~~~~~~~~~~~

   If, in addition to ``PROFILE = TRUE``, you set ``TRACE_PROFILE = TRUE``,
   then the profiler keeps track of when each profiled function is called and
   the ``bl_prof`` directory will include the function call stack. This is
   especially useful when core functions, such as :cpp:`FillBoundary` can be
   called from many different regions of the code. Part of the trace profiling
   is the ability to set regions in the code which can be analyzed for
   profiling information independently from other regions.

Communication Profiling
~~~~~~~~~~~~~~~~~~~~~~~

  If, in addition to ``PROFILE = TRUE``, you set ``COMM_PROFILE = TRUE``, then
  the ``bl_prof`` directory will contain additional information about MPI
  communication (point-to-point timings, data volume, barrier/reduction times,
  etc.). ``TRACE_PROFILE = TRUE`` and ``COMM_PROFILE = TRUE`` can be set
  together.

The AMReX-specific profiling tools are currently under development and this
documentation will reflect the latest status in the development branch.

Instrumenting C++ Code
======================

.. highlight:: c++

You must at least instrument main(), i.e

::

    int main(...)
    {
      amrex::Initialize(argc,argv);
      BL_PROFILE_VAR("main()",pmain);

      ...

      BL_PROFILE_VAR_STOP(pmain);
      amrex::Finalize();
    }

You can then instrument any of your functions

::

    void YourClass::YourFunction()
    {
      BL_PROFILE_VAR("YourClass::YourFunction()",object_name);  // this name can be any string

      // your function code
    }

Note that you do not need to put BL_PROFILE_VAR_STOP because the profiler will
go out of scope at the end of the function.

For other timers within an already instrumented function, add:

::

          BL_PROFILE_VAR("Flaten::FORT_FLATENX()", anyname);  // add this before
            FORT_FLATENX(arg1, arg2);
          BL_PROFILE_VAR_STOP(anyname);   // add this after, using the same name

if you want to use the same name within the same scope, you can use:

::

          BL_PROFILE_VAR("MyFuncs()", myfuncs);  // the first one
            MyFunc_0(arg);
          BL_PROFILE_VAR_STOP(myfuncs);
          ...
          BL_PROFILE_VAR_START(myfuncs);
            MyFunc_1(arg);
          BL_PROFILE_VAR_STOP(myfuncs);

or create a profiling variable without starting, then start/stop:

::

          BL_PROFILE_VAR_NS("MyFuncs()", myfuncs);  // dont start the timer
          ...
          BL_PROFILE_VAR_START(myfuncs);
            MyFunc_0(arg);
          BL_PROFILE_VAR_STOP(myfuncs);
          ...
          BL_PROFILE_VAR_START(myfuncs);
            MyFunc_1(arg);
          BL_PROFILE_VAR_STOP(myfuncs);

Instrumenting Fortran90 Code
============================

When using the full profiler, Fortran90 functions can also be instrumented
with the following calls:

.. highlight:: fortran

::

    call bl_proffortfuncstart("my_function")
    ...
    call bl_proffortfuncstop("my_function")

Note that the start and stop calls must be matched and the profiling output
will warn of any :fortran:`bl_proffortfuncstart` calls that were not stopped
with :fortran:`bl_proffortfuncstop` calls (in debug mode only). You will need
to add :fortran:`bl_proffortfuncstop` before any returns and at the end of the
function or at the point in the function you want to stop profiling.

For functions with a high number of calls, there is a lighter-weight interface:

::

     call bl_proffortfuncstart_int(n)
     ...
     call bl_proffortfuncstop_int(n)

where ``n`` is an integer in the range ``[1,mFortProfsIntMaxFuncs]``.
``mFortProfsIntMaxFuncs`` is currently set to 32.  The profiled
function will be named ``FORTFUNC_n`` in the profiler output,
unless you rename it with ``BL_PROFILE_CHANGE_FORT_INT_NAME(fname, int)``
where ``fname`` is a std::string and ``int`` is the integer ``n``
in the ``bl_proffortfuncstart_int/bl_proffortfuncstop_int`` calls.
``BL_PROFILE_CHANGE_FORT_INT_NAME`` should be called in ``main()``.

Be aware: Fortran functions cannot be profiled when using the Tiny Profiler.
You will need to turn on the full profiler to recieve the results from
fortran instrumentation.

.. _sec:sample:tiny:

Sample Output From Tiny Profile
===============================

Sample output using ``TINY_PROFILE = TRUE`` can look like the following:

.. highlight:: console

::


    TinyProfiler total time across processes [min...avg...max]: 1.765...1.765...1.765
    ---------------------------------------------------------------------------------
    Name                          NCalls   Excl. Min   Excl. Avg   Excl. Max   Max  %
    ---------------------------------------------------------------------------------
    mfix_level::EvolveFluid       1        1.602       1.668       1.691       95.83%
    FabArray::FillBoundary()      11081    0.02195     0.03336     0.06617      3.75%
    FabArrayBase::getFB()         22162    0.02031     0.02147     0.02275      1.29%
    PC<...>::WriteAsciiFile()     1        0.00292     0.004072    0.004551     0.26%


    ---------------------------------------------------------------------------------
    Name                          NCalls   Incl. Min   Incl. Avg  Incl. Max    Max  %
    ---------------------------------------------------------------------------------
    mfix_level::Evolve()          1        1.69        1.723      1.734        98.23%
    mfix_level::EvolveFluid       1        1.69        1.723      1.734        98.23%
    FabArray::FillBoundary()      11081    0.04236     0.05485    0.08826       5.00%
    FabArrayBase::getFB()         22162    0.02031     0.02149    0.02275       1.29%

AMRProfParser
=============

:cpp:`AMRProfParser` is a tool for processing and analyzing the ``bl_prof``
database. It is a command line application that can create performance
summaries, plotfiles showing point to point communication and timelines, HTML
call trees, communication call statistics, function timing graphs, and other
data products. The parser's data services functionality can be called from an
interactive environment such as Amrvis, from a sidecar for dynamic performance
optimization, and from other utilities such as the command line version of the
parser itself. It has been integrated into Amrvis for visual interpretation of
the data allowing Amrvis to open the bl_prof database like a plotfile but with
interfaces appropriate to profiling data. AMRProfParser and Amrvis can be run
in parallel both interactively and in batch mode.

