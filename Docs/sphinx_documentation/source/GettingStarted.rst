.. role:: cpp(code)
   :language: c++


Downloading the Code
====================

The source code is available at https://github.com/AMReX-Codes/amrex. The
GitHub repo is our central repo for development. The development branch
includes the latest state of the code, and it is tagged as a release
on a monthly basis with version number YY.MM (e.g., 17.04). The MM part of the
version is incremented every month, and the YY part every year.  Bug fix
releases are tagged with YY.MM.patch (e.g., 17.04.1).

Each release on GitHub carries a source tarball, ``amrex-YY.MM.tar.gz``,
together with a signed build-provenance attestation,
``amrex-YY.MM.intoto.jsonl``, that records the commit it was made from. To
check a downloaded tarball with the GitHub CLI, run

.. code-block:: bash

   gh attestation verify amrex-YY.MM.tar.gz --repo AMReX-Codes/amrex

The release also carries ``amrex-YY.MM.spdx.json``, an SBOM that lists the
GitHub Actions used by AMReX's CI rather than the library's build
dependencies.

AMReX can also be obtained using Spack (https://spack.io/).  Assuming
you have Spack installed, simply type, ``spack install amrex``. For more
information see the :ref:`sec:build:spack` section in Building AMReX.

Example: Hello World
====================

The source code of this example is at ``amrex-tutorials/Base/HelloWorld/``
and is also shown below.

.. highlight:: c++

::

    #include <AMReX.H>
    #include <AMReX_Gpu.H>
    #include <AMReX_OpenMP.H>
    #include <AMReX_ParallelDescriptor.H>
    #include <AMReX_Print.H>

    int main (int argc, char* argv[])
    {
        // Initialize MPI and the GPU device, if enabled.
        amrex::Initialize(argc, argv);
        {
            // amrex::Print only prints on the I/O process (rank 0).
            amrex::Print() << "Hello world from AMReX version "
                           << amrex::Version() << "\n";

            amrex::Print() << "  MPI ranks:      " << amrex::ParallelDescriptor::NProcs() << "\n"
                           << "  OpenMP threads: " << amrex::OpenMP::get_max_threads() << "\n"
                           << "  GPU backend:    "
    #if defined(AMREX_USE_CUDA)
                           << "CUDA"
    #elif defined(AMREX_USE_HIP)
                           << "HIP"
    #elif defined(AMREX_USE_SYCL)
                           << "SYCL"
    #else
                           << "none"
    #endif
                           << "\n";

            // amrex::AllPrint prints on every process.
            amrex::AllPrint() << "  Rank " << amrex::ParallelDescriptor::MyProc()
                              << " says hello\n";
        }
        amrex::Finalize();
    }

Usually the first and last statements of the :cpp:`int main(...)` function
of every program should be calling :cpp:`amrex::Initialize` and
:cpp:`amrex::Finalize`, respectively. They set up and shut down MPI and the
GPU device, if enabled. In between, :cpp:`amrex::Print` prints a string that
includes the AMReX version returned by the :cpp:`amrex::Version` function,
the number of MPI processes and OpenMP threads, and the GPU backend.
:cpp:`amrex::Print` prints only on the I/O process, whereas
:cpp:`amrex::AllPrint` prints on every process. AMReX objects must be
destroyed before :cpp:`amrex::Finalize` is called, hence the extra scope.
Note that the name of all AMReX header files starts with ``AMReX_`` (or just
AMReX in the case of AMReX.H). All AMReX C++ functions are in the
:cpp:`amrex` namespace.

Building the Code
-----------------

You build the code in the ``amrex-tutorials/Base/HelloWorld/`` directory.
The ``GNUmakefile`` there expects AMReX at ``../../../amrex``, i.e., cloned
next to ``amrex-tutorials``; otherwise set the ``AMREX_HOME`` variable to
the path of your AMReX source tree. Typing ``make USE_MPI=FALSE`` will
start the compilation process and result in an executable named
``main3d.gnu.ex``. The name shows that the GNU compiler is used and that the
executable is built for 3D. Although this simple example code is dimension
independent, dimensionality does matter for all non-trivial examples. The
build process can be adjusted by modifying the
``amrex-tutorials/Base/HelloWorld/GNUmakefile`` file. More details on how to
build AMReX can be found in :ref:`Chap:BuildingAMReX`. The tutorials can
also be built with CMake; see `Tutorials Getting Started`_ in the tutorials
documentation.

.. _`Tutorials Getting Started`: https://amrex-codes.github.io/amrex/tutorials_html/GettingStarted.html

Running the Code
----------------

The example code can be run as follows,

.. highlight:: console

::

      ./main3d.gnu.ex

The result may look like:

.. highlight:: console

::

      Initializing AMReX (26.11-6-gb33ab585fb22)...
      AMReX (26.11-6-gb33ab585fb22) initialized
      Hello world from AMReX version 26.11-6-gb33ab585fb22
        MPI ranks:      1
        OpenMP threads: 1
        GPU backend:    none
        Rank 0 says hello
      AMReX (26.11-6-gb33ab585fb22) finalized

The version string means the current commit b33ab585fb22 (note that the first
letter g in gb33a.. is not part of the hash) is based on 26.11 with 6
additional commits. If the AMReX working tree is dirty (i.e., there are
uncommitted changes), ``-dirty`` is appended.

In the GNUmakefile there are compilation options for DEBUG mode (less optimized
code with more error checking), dimensionality, compiler type, and flags to
enable MPI and/or OpenMP parallelism.  If there are multiple instances of a
parameter, the last instance takes precedence.

Parallelization
---------------

Now let's build with MPI by typing ``make`` (``USE_MPI=TRUE`` is the default
in the GNUmakefile). This should make an executable named
``main3d.gnu.MPI.ex``. Note MPI in the file name. You can then run,

.. highlight:: console

::

      mpiexec -n 4 ./main3d.gnu.MPI.ex amrex.v=1

The result may look like:

.. highlight:: console

::

      Initializing AMReX (26.11-6-gb33ab585fb22)...
      MPI initialized with 4 MPI processes
      MPI initialized with thread support level 0
      AMReX (26.11-6-gb33ab585fb22) initialized
      Hello world from AMReX version 26.11-6-gb33ab585fb22
        MPI ranks:      4
        OpenMP threads: 1
        GPU backend:    none
        Rank 0 says hello
        Rank 1 says hello
        Rank 2 says hello
        Rank 3 says hello
      AMReX (26.11-6-gb33ab585fb22) finalized

The lines printed by :cpp:`amrex::AllPrint` may appear in any order.
If the compilation fails, you are referred to :ref:`Chap:BuildingAMReX` for
more details on how to configure the build system.  The *optional* command line
argument ``amrex.v=1`` sets the AMReX verbosity level
to 1 to print the number of MPI processes used.  The default verbosity
level is 1, and you can pass ``amrex.v=0`` to turn it off.
More details on how runtime parameters are handled can be found in
section :ref:`sec:basics:parmparse`.

If you want to build with OpenMP, type ``make USE_OMP=TRUE USE_MPI=FALSE``.
This should make an executable named ``main3d.gnu.OMP.ex``. Note OMP in the
file name. Make sure the ``OMP_NUM_THREADS`` environment variable is set on
your system. You can then run,

.. highlight:: console

::

      OMP_NUM_THREADS=4 ./main3d.gnu.OMP.ex

The result may look like:

.. highlight:: console

::

      Initializing AMReX (26.11-6-gb33ab585fb22)...
      OMP initialized with 4 OMP threads
      AMReX (26.11-6-gb33ab585fb22) initialized
      Hello world from AMReX version 26.11-6-gb33ab585fb22
        MPI ranks:      1
        OpenMP threads: 4
        GPU backend:    none
        Rank 0 says hello
      AMReX (26.11-6-gb33ab585fb22) finalized

Note that you can build with both ``USE_MPI=TRUE`` and ``USE_OMP=TRUE``.  You
can then run,

.. highlight:: console

::

      OMP_NUM_THREADS=4 mpiexec -n 2 ./main3d.gnu.MPI.OMP.ex

The result may look like:

.. highlight:: console

::

      Initializing AMReX (26.11-6-gb33ab585fb22)...
      MPI initialized with 2 MPI processes
      MPI initialized with thread support level 0
      OMP initialized with 4 OMP threads
      AMReX (26.11-6-gb33ab585fb22) initialized
      Hello world from AMReX version 26.11-6-gb33ab585fb22
        MPI ranks:      2
        OpenMP threads: 4
        GPU backend:    none
        Rank 0 says hello
        Rank 1 says hello
      AMReX (26.11-6-gb33ab585fb22) finalized

.. _sec:heat equation:

Example: Heat Equation Solver
=============================

We now look at a more complicated example at
``amrex-tutorials/Base/HeatEquation`` and show how simulation results
can be visualized. This example solves the heat equation,

.. math:: \frac{\partial\phi}{\partial t} = D\nabla^2\phi

using forward Euler temporal integration on a periodic domain, with the
standard 5-point (in 2D) or 7-point (in 3D) stencil for the Laplacian, e.g.,
in 2D,

.. math::

   \phi_{i,\,j}^{n+1} = \phi_{i,\,j}^n
   + D\Delta t\left(\frac{\phi_{i+1,\,j}^n-2\phi_{i,\,j}^n+\phi_{i-1,\,j}^n}{\Delta x^2}
   + \frac{\phi_{i,\,j+1}^n-2\phi_{i,\,j}^n+\phi_{i,\,j-1}^n}{\Delta y^2}\right).

The implementation details of the code are discussed in the `Heat Equation`_
tutorial. For now let's just build and run the code, and visualize the
results.

.. _`Heat Equation`: https://amrex-codes.github.io/amrex/tutorials_html/Base.html#base-heatequation

Building and Running the Code
-----------------------------

To build a 2D executable, go to
``amrex-tutorials/Base/HeatEquation`` and type ``make DIM=2``. This
will generate an executable named ``main2d.gnu.MPI.ex``. To run it, type,

.. highlight:: console

::

      mpiexec -n 2 ./main2d.gnu.MPI.ex inputs plot_int=10

Note that the command takes the file ``inputs``, followed by a runtime
parameter that overrides the value in the file. The calculation solves the
heat equation in 2D on a domain with :math:`32 \times 32` cells until
:math:`t = 0.01` and makes a plotfile every :math:`10` steps.  When the run
finishes, you will have a number of plotfiles, ``plt00000, plt00010,`` etc,
in the directory where you are running, and the code prints the error
against the exact solution.  You can control runtime parameters such as the
number of cells, the stop time and how often to write plotfiles by setting
them in ``inputs`` or on the command line.

Visualization
=============

There are several visualization tools that can be used for AMReX plotfiles.
One standard tool used within the AMReX community is Amrvis, a package
developed and supported by CCSE that is designed specifically for highly
efficient visualization of block-structured hierarchical AMR data.  (Amrvis can
also be used to visualize performance data; see the :ref:`Chap:AMRex-based
Profiling Tools` chapter for further details.) Plotfiles can also be viewed
using the VisIt, ParaView, and yt packages.  Particle data can be viewed using
ParaView.  Refer to Chapter on :ref:`Chap:Visualization` for how to use each of
these tools.

Tutorials
=========

Users new to AMReX may be interested in the `AMReX Tutorials`_, small
stand-alone codes organized by topic, from introductory examples of the
basic data structures to complete codes meant as starting points for new
applications. All but the Python tutorials build with CMake or GNU Make,
and each is documented together with its source code.

.. _`AMReX Tutorials`: https://amrex-codes.github.io/amrex/tutorials_html/
