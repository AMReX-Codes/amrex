.. role:: cpp(code)


Frequently Asked Questions
==========================


**Q.** Why am I getting a segmentation fault after my code runs?

**A.** Do you have :cpp:`amrex::initialize(); {` and :cpp:`} amrex::finalize();`
at the beginning and end of your code? For all AMReX commands to function
properly, including to release resources, they need to be contained
between these two curly braces. In the `Initialize and Finalize`_ section,
these commands are discussed further detail.

.. _`Initialize and Finalize` : https://amrex-codes.github.io/amrex/docs_html/Basics.html#initialize-and-finalize

|

**Q.** I want to use a different compiler with GNU Make to compile AMReX. How do I do this?

**A.** In the file ``amrex/Tools/GNUMake/Make.local`` you can specify your own compile
commands by setting the variables ``CXX``, ``CC``, ``FC``, and ``F90``.
An example can be found at `Specifying your own compiler`_ . Additional
customizations are described in the file, ``amrex/Tools/GNUMake/Make.local.template``.
In the same directory, ``amrex/Tools/GNUMake/README.md`` contains detailed
information on compiler commands.

.. _`Specifying your own compiler` : https://amrex-codes.github.io/amrex/docs_html/BuildingAMReX.html#specifying-your-own-compiler

|

**Q.** I'm having trouble compiling my code.

**A.** AMReX developers have found that running the command ``make realclean`` can resolve
many compilation issues.

If you are working in an environment that uses
a module system, please ensure you have the correct modules loaded. Typically, to do this,
type ``module list`` at the command prompt.

|

**Q.** When I profile my code that uses GPUs with ``TINY_PROFILE=TRUE`` or ``PROFILE=TRUE``
my timings are inconsistent.

**A.** Due to the asynchronous nature of GPU execution, data collected by the built-in profilers
may contain some variation. GPU profiling is an area of current development for the AMReX team.

|

**Q.** How do I know I am getting the right answer?

**A.** AMReX provides support for verifying output with several tools. To briefly mention a few:

- The :cpp:`print_state` command can be used to output the data of a single cell.
- :cpp:`VisMF::Write` can be used to write MultiFab data to disk that can be viewed with `Amrvis`_.
- :cpp:`amrex::Print()` is useful for printing
  output when using multiple processes or threads as it prevents messages
  from getting mixed up.
- `fcompare`_ compares two plotfiles and reports absolute and relative error.

Additional tools and discussion on this topic is contained
in the section `Debugging`_.

.. _`Debugging`: https://amrex-codes.github.io/amrex/docs_html/Basics.html#debugging

.. _`Amrvis`: https://amrex-codes.github.io/amrex/docs_html/Visualization.html#sec-amrvis

.. _`fcompare`: https://amrex-codes.github.io/amrex/docs_html/Post_Processing.html#fcompare

|

**Q.** What's the difference between :cpp:`copy` and :cpp:`parallelCopy`?

**A.**

|

**Q.** How do I fill ghost cells?

**A.**

|

**Q.** What's the difference between ``AmrCore`` and ``AmrLevel``? How do
I decide which to use?

**A.**

|

**Q.** For GPU usage, how can I perform explicit host to device and
device to host copies without relying on managed memory?

**A.** Use ``The_Pinned_Arena()`` and ``htod_memcpy()`` or ``dtoh_memcpy()``. See
`Memory Allocation`_ in the AMReX Source Documentation.

.. _`Memory Allocation`: https://amrex-codes.github.io/amrex/docs_html/GPU.html#memory-allocation

|

**Q.** How can you prevent a section of code from running on the GPU for profiling purposes?

**A.**

|

**Q.** How do I generate random numbers with AMReX? Can I set the seed?
Are they thread safe with MPI and OpenMP?

**A.** (Thread safe) Yes, :cpp:`amrex::Random()` is thread safe. When OpenMP is on,
each thread will have its own dedicated Random Number Generator that
is totally independent of the others.
|

**Q.** Is Dirichlet boundary condition data loaded into cell-centered, or
face-centered containers? How is it used in AMReX-based codes like MLMG and the
advection routines in AMReX-Hydro?

**A.** In the cell-centered MLMG solver, the Dirichlet boundary data are stored
in containers that have the information of the location of the data.

|

**Q.** When using embedded boundaries (EB), is :cpp:`flag.isRegular()` the same
as :cpp:`volfrac==1`?

**A.**

|

**Q.** When using embedded boundaries (EB), how far out does
:cpp:`flag.isConnected(ii,jj,kk)` go? How does a cell ``(i,j,k)``
know if a cell ``(i+1,j+1,k+1)`` is "connected" to it?

**A.**

|
|

More Questions
--------------

If your question was not addressed here, you are encouraged to
search and ask for help on the `AMReX GitHub Discussions`_ page.

.. _`AMReX GitHub Discussions`: https://github.com/AMReX-Codes/amrex/discussions


