.. role:: cpp(code)


Frequently Asked Questions
==========================


**Q.** Why am I getting a segmentation fault after my code runs?

**A.** Do you have :cpp:`amrex::initialize(); {` and :cpp:`} amrex::finalize();`
at the beginning and end of your code? For all AMReX commands to function
properly, including to release resources, they need to be contained
between these two curly braces or in a separate function. In the `Initialize
and Finalize`_ section, these commands are discussed further detail.

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

**A.** AMReX developers have found that running the command ``make clean`` can resolve
many compilation issues.

If you are working in an environment that uses
a module system, please ensure you have the correct modules loaded. Typically, to do this,
type ``module list`` at the command prompt.

|

**Q.** When I profile my code that uses GPUs with ``TINY_PROFILE=TRUE`` or ``PROFILE=TRUE``
my timings are inconsistent.

**A.** Due to the asynchronous nature of GPU execution, profilers might only
measure the run time on CPU, if there is no explicit synchronization.  For
``TINY_PROFILE``, one could use :cpp:`ParmParse` parameter
``tiny_profiler.device_synchronize_around_region=1`` to add synchronization.
Note that this may degrade performance.

|

**Q.** How do I know I am getting the right answer?

**A.** AMReX provides support for verifying output with several tools. To briefly mention a few:

- The :cpp:`print_state` function can be used to output the data of a single cell.
- :cpp:`VisMF::Write` can be used to write MultiFab data to disk that can be viewed with `Amrvis`_.
- :cpp:`amrex::Print()` and :cpp:`amrex::AllPrint()` are useful for printing
  output when using multiple processes or threads as it prevents messages
  from getting mixed up.
- `fcompare`_ compares two plotfiles and reports absolute and relative error.

Additional tools and discussion on this topic is contained
in the section `Debugging`_.

.. _`Debugging`: https://amrex-codes.github.io/amrex/docs_html/Basics.html#debugging

.. _`Amrvis`: https://amrex-codes.github.io/amrex/docs_html/Visualization.html#sec-amrvis

.. _`fcompare`: https://amrex-codes.github.io/amrex/docs_html/Post_Processing.html#fcompare

|

**Q.** What's the difference between :cpp:`Copy` and :cpp:`ParallelCopy` for
:cpp:`MultiFab` data?

**A.** :cpp:`MultiFab::Copy` is for two :cpp:`MultiFab`s built with the same
:cpp:`BoxArray` and :cpp:`DistributionMapping`, whereas :cpp:`ParallelCopy`
is for parallel communication of two :cpp:`MultiFab`s with different
:cpp:`BoxArray` and/or :cpp:`DistributionMapping`.

|

**Q.** How do I fill ghost cells?

**A.** See `Ghost Cells`_ in the AMReX Source Documentation.

.. _`Ghost Cells`: https://amrex-codes.github.io/amrex/docs_html/Basics.html#ghost-cells

|

**Q.** What's the difference between ``AmrCore`` and ``AmrLevel``? How do
I decide which to use?

**A.** The :cpp:`AmrLevel` class is an abstract base class that holds data
for a single AMR level.  A vector of :cpp:`AmrLevel` is stored in the
:cpp:`Amr` class, which is derived from :cpp:`AmrCore`.  An application code
can derive from :cpp:`AmrLevel` and override functions.  :cpp:`AmrCore`
contains the meta-data for the AMR hierarchy, but it does not contain any
floating point mesh data.  Instead of using :cpp:`Amr`/:cpp:`AmrLevel`, an
application can also derive from :cpp:`AmrCore`.  If you want flexibility,
you might choose the :cpp:`AmrCore` approach, otherwise the :cpp:`AmrLevel`
approach might be easier because it already has a lot of built-in
capabilities that are common for AMR applications.

|

**Q.** For GPU usage, how can I perform explicit host to device and
device to host copies without relying on managed memory?

**A.** Use ``The_Pinned_Arena()`` (See `Memory Allocation`_ in the AMReX
Source Documentation.) and

.. code-block::

 void htod_memcpy (void* p_d, const void* p_h, const std::size_t sz);
 void dtoh_memcpy (void* p_h, const void* p_d, const std::size_t sz);
 void dtoh_memcpy (FabArray<FAB>& dst, FabArray<FAB> const& src, int scomp, int dcomp, int ncomp);
 void htod_memcpy (FabArray<FAB>& dst, FabArray<FAB> const& src, int scomp, int dcomp, int ncomp);

.. _`Memory Allocation`: https://amrex-codes.github.io/amrex/docs_html/GPU.html#memory-allocation

|

**Q.** How can I prevent a section of code from running on the GPU?

**A.** Use:

.. code-block::

    Gpu::setLaunchRegion(0);
    ... ;
    Gpu::setLaunchRegion(1);

Please note that because much of the execution patterns remain intact with this approach,
it is likely not the ideal way to compare GPU to non-GPU performance. For more information
see `Basic Gpu Debugging`_.

.. _`Basic Gpu Debugging`: GPU.html#basic-gpu-debugging

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
