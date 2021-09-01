.. role:: cpp(code)
   :language: c++

.. role:: fortran(code)
   :language: fortran

Types of Profiling
==================

AMReX's build-in profiling works through objects that start and stop timers
based on user-placed macros or in the objects' constructor and destructor.
The results of these timers are stored in a global list that is consolidated
and printed during finalization, or at a user-defined flush.

Currently, AMReX has two options for build-in profiling:
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

  AMReX_TINY_PROFILE = ON
  AMReX_BASE_PROFILE = OFF

Note that if you set ``PROFILE = TRUE``  (or ``AMReX_BASE_PROFILE =
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

AMReX profiler objects are created and managed through :cpp:`BL_PROF` macros.

.. highlight:: c++

To start with, you must at least instrument main(), i.e.:

::

    int main(...)
    {
      amrex::Initialize(argc,argv);
      BL_PROFILE_VAR("main()",pmain);

      <AMReX code block>

      BL_PROFILE_VAR_STOP(pmain);
      amrex::Finalize();
    }

    // OR

    void main_main()
    {
        BL_PROFILE("main()");

        <AMReX code block>
    }

    int main(...)
    {
        amrex::Initialize(argc,argv);
        main_main();
        amrex::Finalize();
    }

You can then instrument any of your functions, or code blocks. There are four general
profiler macro types available:

1) A scoped timer, :cpp:`BL_PROFILE`:

These timers generate their own object name, so they can't be controlled after defined.
However, they are the cleanest and easiest to work with in many situations. They time from
the point the macro is called until the end of the enclosing scope. This macro is ideal to
time an entire function. For example:

::

    void YourClass::YourFunction()
    {
      BL_PROFILE("YourClass::YourFunction()");   // Timer starts here.

      < Your Function Code Block>

    }    // <------ Timer goes out of scope here, calling stop and returning the function time.

Note that all AMReX timers are scoped and will call "stop" when the corresponding object is destroyed.
This macro is unique because it can _only_ stop when it goes out of scope.

2) A named, scoped timer, :cpp:`BL_PROFILE_VAR`:

In some cases, using scopes to control the timer is non-ideal. In such cases, you can use the
`_VAR_` macros to create a named timer that can be controlled through `_START_` and `_STOP_` macros.
`_VAR_` signifies the macro takes a variable name. For example, to time a function without scoping:

::
          BL_PROFILE_VAR("Flaten::FORT_FLATENX()", anyname);  // Create and start "anyname".
            FORT_FLATENX(arg1, arg2);
          BL_PROFILE_VAR_STOP(anyname);   // Stop the "anyname" timer object.

This can also be used to selectively time with the same scope, for example, to include :cpp:`Func_0`
and :cpp:`Func_2`, but not :cpp:`Func_1`:

::

          BL_PROFILE_VAR("MyFuncs()", myfuncs);  // the first one
            MyFunc_0(args);
          BL_PROFILE_VAR_STOP(myfuncs);

            MyFunc_1(args);

          BL_PROFILE_VAR_START(myfuncs);
            MyFunc_2(arg);
          BL_PROFILE_VAR_STOP(myfuncs);

Remember, these are still scoped. So, the scoped timer example can be reproduced exactly with named
timers by just using the :cpp:`_VAR` macro:

::

    void YourClass::YourFunction()
    {
      BL_PROFILE_VAR("YourClass::YourFunction()",  pmain);   // Timer starts here.

      < Your Function Code Block>

    }    // <------ Timer goes out of scope here correctly, without a STOP call.



3) A named, scoped timer that doesn't auto-start, :cpp:`BL_PROFILE_VAR_NS`:

Sometimes, a complicated scoping may mean the profiling object needs to be defined before it's
started. To create a named AMReX timer that doesn't start automatically, use the `_NS_` macros.
(NS stands for "no start").

For example, this implementation times :cpp:`MyFunc0` and :cpp:`MyFunc1` but not any of the
"Additional Code" blocks:

::

          {
              BL_PROFILE_VAR_NS("MyFuncs()", myfuncs);  // dont start the timer

              <Additional Code A>

              {
                 BL_PROFILE_VAR_START(myfuncs);
                   MyFunc_0(arg);
                 BL_PROFILE_VAR_STOP(myfuncs);
              }

              <Additional Code B>

              {
                 BL_PROFILE_VAR_START(myfuncs);
                   MyFunc_1(arg);
                 BL_PROFILE_VAR_STOP(myfuncs);

                 <Additional Code C>
              }
          }

Note: The `_NS_` macro must, by necessity, also be a `_VAR_` macro, otherwise, you would never be
able to turn the timer on!

4) Designate a sub-region to profile, :cpp:`BL_PROFILE_REGION`:

Often, it's helpful to look at a subset of timers separately from the complete profile. For
example, look at the timing of a specific timestep or isolate everything inside the "Chemistry"
part of the code. This can be accomplished by designating profile regions. All timers within a
named region will be included both in the full analysis, as well as a separate sub-analysis for
further analysis.

Regions are meant to be large, contiguous blocks of code and should be used sparingly and purposefully
to produce useful profiling report. As such, the possible region options are purposefully limited.
When using the TinyProfiler, the only available region macro is the scoped macro. To create a region
that includes in the `MyFuncs` code block, including timers in the "Additional Code" regions:

::

          {
              BL_PROFILE_REGION("MyFuncs");

              <Additional Code A>

              {
                 BL_PROFILE("MyFunc0");

                 MyFunc_0(arg);
              }

              <Additional Code B>

              {
                 BL_PROFILE("MyFunc1");

                 MyFunc_1(arg);
                 <Additional Code C>
              }
          }

If using the Full Profiler, named region objects are also available. These use a slightly modified
`_VAR_`, `_START_` and `_STOP_` formatting. To include all timers in a region, except the timers
inside :cpp:`MyFunc0` and :cpp:`MyFunc1`:

::
          {
              BL_PROFILE_REGION("MyFuncs()", myfuncs);
              <Code Block A>
              BL_PROFILE_REGION_STOP(myfuncs);

              {
                   MyFunc_0(arg);
              }

              BL_PROFILE_REGION_START(myfuncs);
              <Code Block B>
              BL_PROFILE_REGION_STOP(myfuncs);

              {
                   MyFunc_1(arg);

                 BL_PROFILE_REGION_START(myfuncs);
                 <Code Block C>
                 BL_PROFILE_REGION_STOP(myfuncs);
              }
          }

Currently, these are only available in the FullProfiler to control the size of the TinyProfiler output
and allow instrumentation strategies without needing to re-code each time. You can use the scoped
regions in a few select places that will be output in both profiler modes, while named regions
will only be isolated during a FullProfiler run.


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
You will need to turn on the full profiler to receive the results from
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
