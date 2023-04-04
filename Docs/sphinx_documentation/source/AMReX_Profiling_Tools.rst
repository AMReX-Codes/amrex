.. role:: cpp(code)
   :language: c++

.. role:: fortran(code)
   :language: fortran

Types of Profiling
==================

AMReX's built-in profiling works through objects that start and stop timers
based on user-placed macros or an object's constructor and destructor.
The results from these timers are stored in a global list that is consolidated
and printed during finalization, or at a user-defined flush point.

Currently, AMReX has two options for built-in profiling:
:ref:`sec:tiny:profiling` and :ref:`sec:full:profiling`.

.. _sec:tiny:profiling:

Tiny Profiling
----------------------

To enable "Tiny Profiling" with GNU Make edit the options in the file ``GNUMakefile``
to show,

::

  TINY_PROFILE = TRUE
  PROFILE      = FALSE

If building with CMake, set the following CMake flags,

::

  AMReX_TINY_PROFILE = ON
  AMReX_BASE_PROFILE = OFF

.. note::
   If you set ``PROFILE = TRUE`` (or ``AMReX_BASE_PROFILE =
   ON``) to enable full profiling then this will override the ``TINY_PROFILE`` flag
   and tiny profiling will be disabled.

Output
~~~~~~

At the end of a run, a summary of exclusive and inclusive function times will
be written to ``stdout``. This output includes the minimum and maximum (over
processes) time spent in each routine as well as the average and the maximum
percentage of total run time. See the sample output below.

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


The tiny profiler automatically writes the results to ``stdout`` at the end of your
code, when ``amrex::Finalize();`` is reached. However, you may want to write
partial profiling results to ensure your information is saved when you may fail
to converge or if you expect to run out of allocated time. Partial results can
be written at user-defined points in the code by inserting the line:

::

  BL_PROFILE_TINY_FLUSH();

Any timers that have not reached their ``BL_PROFILE_VAR_STOP`` call or exited
their scope and deconstructed will not be included in these partial outputs.
(e.g., a properly instrumented ``main()`` should show a time of zero in all
partial outputs.) Therefore, it is recommended to place these flush calls in
easily identifiable regions of your code and outside of as many profiling
timers as possible, such as immediately before or after writing a checkpoint.

Also, since flush calls will print multiple, similar looking outputs to ``stdout``,
it is also recommended to wrap any ``BL_PROFILE_TINY_FLUSH();`` calls in
informative ``amrex::Print()`` lines to ensure accurate identification of each
set of timers.

.. _sec:full:profiling:

Full Profiling
--------------

If you set ``PROFILE = TRUE`` then a ``bl_prof`` directory will be written that
contains detailed per-task timings for each processor.  This will be written in
``nfiles`` files (where ``nfiles`` is specified by the user). The information
in the directory can be analyzed by the :ref:`sec:amrprofparse` tool
within :ref:`sec:amrvis`. In addition, an
exclusive-only set of function timings will be written to ``stdout``.

Trace Profiling
~~~~~~~~~~~~~~~

   If you set ``TRACE_PROFILE = TRUE`` in addition to ``PROFILE = TRUE``,
   then the profiler keeps track of when each profiled function is called and
   the ``bl_prof`` directory will include the function call stack. This is
   especially useful when core functions, such as :cpp:`FillBoundary` can be
   called from many different regions of the code. Using trace profiling
   allows one to specify regions in the code that can be analyzed for
   profiling information independently from other regions.

Communication Profiling
~~~~~~~~~~~~~~~~~~~~~~~

  If you set ``COMM_PROFILE = TRUE`` in addition to ``PROFILE = TRUE``, then
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

To start, you must at least instrument main(), i.e.:

::

    int main(...)
    {
      amrex::Initialize(argc,argv);
      BL_PROFILE_VAR("main()",pmain);

      <AMReX code block>

      BL_PROFILE_VAR_STOP(pmain);
      amrex::Finalize();
    }

Or:

::

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
----------------------------------------------------------------------------------------

These timers generate their own object names, so they can't be controlled after being defined.
However, they are the cleanest and easiest to work with in many situations. They time from
the point where the macro is called until the end of the enclosing scope. This macro is ideal
for timing an entire function. For example:

::

    void YourClass::YourFunction()
    {
      BL_PROFILE("YourClass::YourFunction()");   // Timer starts here.

      < Your Function Code Block>

    }    // <------ Timer goes out of scope here, calling stop and returning the function time.

Note that all AMReX timers are scoped and will call "stop" when the corresponding object is destroyed.
This macro is unique because it can *only* stop when it goes out of scope.

2) A named, scoped timer, :cpp:`BL_PROFILE_VAR`:
----------------------------------------------------------------------------------------

In some cases, using scopes to control a timer is not ideal. In such cases, you can use the
``_VAR_`` macros to create a named timer that can be controlled through ``_START_`` and ``_STOP_`` macros.
``_VAR_`` signifies that the macro takes a variable name. For example, to time a function without scoping:

::

          BL_PROFILE_VAR("Flaten::FORT_FLATENX()", anyname);  // Create and start "anyname".
            FORT_FLATENX(arg1, arg2);
          BL_PROFILE_VAR_STOP(anyname);   // Stop the "anyname" timer object.

This can also be used to selectively time with the same scope. For example, to include :cpp:`Func_0`
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
----------------------------------------------------------------------------------------

Sometimes, a complicated scoping may mean the profiling object needs to be defined before it's
started. To create a named AMReX timer that doesn't start automatically, use the ``_NS_`` macros.
("NS" stands for "no start"). For example, this implementation times :cpp:`MyFunc0`
and :cpp:`MyFunc1` but not any of the
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

.. Note::
   The ``_NS_`` macro must, by necessity, also be a ``_VAR_`` macro. Otherwise, you would never be
   able to turn the timer on!

4) Designate a sub-region to profile, :cpp:`BL_PROFILE_REGION`:
----------------------------------------------------------------------------------------

Often, it's helpful to look at a subset of timers separately from the complete profile. For
example, you may want to view the timing of a specific time step or isolate everything inside the "Chemistry"
part of the code. This can be accomplished by designating profile regions. All timers within a
named region will be included both in the full analysis, as well as in a separate sub-analysis.

Regions are meant to be large contiguous blocks of code, and should be used sparingly and purposefully
to produce useful profiling reports. As such, the possible region options are purposefully limited.

Scoped Regions
~~~~~~~~~~~~~~

When using the Tiny Profiler, the only available region macro is the scoped macro. To create a region
that profiles the `MyFuncs` code block, including all timers in the "Additional Code" regions, add
macros in the following way:

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

The ``MyFuncs`` region appears in the Tiny Profiler output as an additional table.
The following output example, mimics the above code. In it, the region is
indicated by ``REG::MyFuncs``.

.. code-block:: console

    BEGIN REGION MyFuncs

    -------------------------------------------------------------
    Name          NCalls  Excl. Min  Excl. Avg  Excl. Max   Max %
    -------------------------------------------------------------
    MyFunc0         1000      4.402      4.402      4.402  14.19%
    MyFunc1         1000       4.39       4.39       4.39  14.15%
    REG::MyFuncs    1000     0.0168     0.0168     0.0168   0.05%
    -------------------------------------------------------------

    -------------------------------------------------------------
    Name          NCalls  Incl. Min  Incl. Avg  Incl. Max   Max %
    -------------------------------------------------------------
    REG::MyFuncs    1000      8.809      8.809      8.809  28.39%
    MyFunc0         1000      4.402      4.402      4.402  14.19%
    MyFunc1         1000       4.39       4.39       4.39  14.15%
    -------------------------------------------------------------

    END REGION MyFuncs


Named Regions
~~~~~~~~~~~~~~

If using the Full Profiler, named region objects are also available.
Named regions allow control of start and stop points without relying on scope.
These macros use slightly modified
``_VAR_``, ``_START_`` and ``_STOP_`` formatting. The first
argument is the name, followed by the profile variable. Names
for each section can differ, but because the profiler variable will be used
to group the sections into a region, it must be the same.
Consider the following example:

::

  {
      BL_PROFILE_REGION_VAR("RegionAC",reg_ac);
      <Code Block A>
      BL_PROFILE_REGION_VAR_STOP("RegionAC", reg_ac);

      {

         MyFunc_0(arg);
      }

      BL_PROFILE_REGION_VAR("RegionB", reg_b)
      <Code Block B>
      BL_PROFILE_REGION_VAR_STOP("RegionB", reg_b);

      {

         MyFunc_1(arg);

         BL_PROFILE_REGION_VAR_START("SecondRegionAC", reg_ac);
         <Code Block C>
         BL_PROFILE_REGION_VAR_STOP("SecondRegionAC", reg_ac);
      }
  }


Here, :cpp:`<Code Block A>` and :cpp:`<Code Block C>` are
grouped into one region labeled "RegionAC" for profiling. :cpp:`<Code Block B>`
is isolated in its own group.
Any timers inside :cpp:`MyFunc_0` and :cpp:`MyFunc_1` are not included in the
region groupings.

Instrumenting Fortran90 Code
============================

When using the full profiler, Fortran90 functions can also be instrumented
with the following calls:

.. highlight:: fortran

::

    call bl_proffortfuncstart("my_function")
    ...
    call bl_proffortfuncstop("my_function")

Note that the start and stop calls must be matched before leaving the
scope of the corresponding start. Moreover, it is necessary to take into
account all possible code paths. Therefore, you may need to add :fortran:`bl_proffortfuncstop`
in multiple locations, such as before any returns, at the end of the function
and at the point in the function where you want to stop profiling. The profiling
output will only warn of any :fortran:`bl_proffortfuncstart` calls that were not stopped with
:fortran:`bl_proffortfuncstop` calls when in debug mode.

For functions with a high number of calls, there is a lighter-weight interface,

::

     call bl_proffortfuncstart_int(n)
     ...
     call bl_proffortfuncstop_int(n)

where ``n`` is an integer in the range ``[1,mFortProfsIntMaxFuncs]``.
``mFortProfsIntMaxFuncs`` is currently set to 32.  The profiled
function will be named ``FORTFUNC_n`` in the profiler output,
unless you rename it with ``BL_PROFILE_CHANGE_FORT_INT_NAME(fname, int)``
where ``fname`` is a ``std::string`` and ``int`` is the integer ``n``
in the ``bl_proffortfuncstart_int/bl_proffortfuncstop_int`` calls.
``BL_PROFILE_CHANGE_FORT_INT_NAME`` should be called in ``main()``.

.. warning::
   Fortran functions cannot be profiled when using the Tiny Profiler.
   You will need to turn on the Full Profiler to receive the results from
   fortran instrumentation.

.. _sec:profopts:

Profiling Options
=================

AMReX's communication algorithms are often regions of code that increase in wall clock time
when the application is load imbalanced, due to the MPI_Wait calls in these functions.
To better understand if this is occuring and by how much, you can turn on an AMReX timed
synchronization with the runtime variable: ``amrex.use_profiler_syncs=1`` This adds named timers
beginning with ``SyncBeforeComms`` immediately prior to the start of the FillBoundary,
ParallelCopy and particle Redistribute functions, isolating any prior load imbalance to that timer
before beginning the comm operation.

This is a diagnostic tool and may slow your code down, so it is not recommended to turn this
on for production runs.

.. note::
  Note: the ``SyncBeforeComms`` timer is not equal to your load imbalance. It only captures imbalance
  between the comm functions and the previous sync point; there may be other load imbalances
  captured elsewhere. Also, the timer reports in terms of MPI rank, so if the most imbalanced
  rank changes throughout the simulation, the timer will be an underestimation.

  The effect on the communication timers may be more helpful: they will show the time to complete
  communications if there was no load imbalance. This means the difference between a case
  with and without this profiler sync may be a more useful metric for analysis.

.. _sec:amrprofparse:

AMRProfParser
=============

:cpp:`AMRProfParser` is a tool for processing and analyzing the ``bl_prof``
database. It is a command line application that can create performance
summaries, plotfiles showing point-to-point communication and timelines, HTML
call trees, communication call statistics, function timing graphs, and other
data products. The parser's data services functionality can be called from an
interactive environment such as :ref:`sec:amrvis`, from a sidecar for dynamic performance
optimization, and from other utilities such as the command line version of the
parser itself. It has been integrated into Amrvis for visual interpretation of
the data allowing Amrvis to open the ``bl_prof`` database like a plotfile but with
interfaces appropriate to profiling data. AMRProfParser and Amrvis can be run
in parallel both interactively and in batch mode.
