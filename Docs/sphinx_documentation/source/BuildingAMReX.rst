.. role:: cpp(code)
   :language: c++


.. _sec:build:make:

Building with GNU Make
======================

In this build approach, you write your own make files defining a number of
variables and rules. Then you invoke  ``make`` to start the building process.
This will result in an executable upon successful completion. The temporary
files generated in the building process are stored in a temporary directory
named  ``tmp_build_dir``.

Dissecting a Simple Make File
-----------------------------

An example of building with GNU Make can be found in
``amrex/Tutorials/Basic/HelloWorld_C``.  :numref:`tab:makevars` below shows a
list of important variables.

.. raw:: latex

   \begin{center}

.. _tab:makevars:

.. table:: Important make variables

   +------------+-------------------------------------+-------------+
   | Variable   | Value                               | Default     |
   +============+=====================================+=============+
   | AMREX_HOME | Path to amrex                       | environment |
   +------------+-------------------------------------+-------------+
   | COMP       | gnu, cray, ibm, intel, llvm, or pgi | none        |
   +------------+-------------------------------------+-------------+
   | DEBUG      | TRUE or FALSE                       | TRUE        |
   +------------+-------------------------------------+-------------+
   | DIM        | 1 or 2 or 3                         | none        |
   +------------+-------------------------------------+-------------+
   | USE_MPI    | TRUE or FALSE                       | FALSE       |
   +------------+-------------------------------------+-------------+
   | USE_OMP    | TRUE or FALSE                       | FALSE       |
   +------------+-------------------------------------+-------------+

.. raw:: latex

   \end{center}

At the beginning of ``amrex/Tutorials/Basic/HelloWorld_C/GNUmakefile``,
``AMREX_HOME`` is set to the path to the top directory of AMReX.  Note that in
the example :cpp:`?=` is a conditional variable assignment operator that only
has an effect if ``AMREX_HOME`` has not been defined (including in the
environment). One can also set ``AMREX_HOME`` as an environment variable. For
example in bash, one can set

.. highlight:: bash

::

    export AMREX_HOME=/path/to/amrex

alternatively, in tcsh one can set

.. highlight:: bash

::

    setenv AMREX_HOME /path/to/amrex

One must set the ``COMP`` variable to choose a compiler. Currently the list of
supported compilers includes gnu, cray, ibm, intel, llvm, and pgi. One must
also set the ``DIM`` variable to either 1, 2, or 3, depending on the dimensionality
of the problem.

Variables ``DEBUG``, ``USE_MPI`` and ``USE_OMP`` are optional with default set
to TRUE, FALSE and FALSE, respectively.  The meaning of these variables should
be obvious.  When ``DEBUG = TRUE``, aggressive compiler optimization flags are
turned off and assertions in  source code are turned on. For production runs,
``DEBUG`` should be set to FALSE.

After defining these make variables, a number of files, ``Make.defs,
Make.package`` and ``Make.rules``, are included in the GNUmakefile. AMReX-based
applications do not need to include all directories in AMReX; an application
which does not use particles, for example, does not need to include files from
the Particle directory in its build.  In this simple example, we only need to
include ``$(AMREX_HOME)/Src/Base/Make.package``. An application code also has
its own Make.package file (e.g., ``./Make.package`` in this example) to append
source files to the build system using operator ``+=``. Variables for various
source files are shown below.

    CEXE_sources
        C++ source files. Note that C++ source files are assumed to have a .cpp
        extension.

    CEXE_headers
        C++ headers with .h or .H extension.

    cEXE_sources
        C source files with .c extension.

    cEXE_headers
        C headers with .h extension.

    f90EXE_sources
        Free format Fortran source with .f90 extension.

    F90EXE_sources
        Free format Fortran source with .F90 extension.  Note that these
        Fortran files will go through preprocessing.

In this simple example, the extra source file, ``main.cpp`` is in the current
directory that is already in the build system's search path. If this example
has files in a subdirectory (e.g., ``mysrcdir``), you will then need to add the
following to ``Make.package``.

::

        VPATH_LOCATIONS += mysrcdir
        INCLUDE_LOCATIONS += mysrcdir

Here ``VPATH_LOCATIONS`` and ``INCLUDE_LOCATIONS`` are the search path for
source and header files, respectively.

Tweaking the Make System
------------------------

The GNU Make build system is located at ``amrex/Tools/GNUMake``.  You can read
``README.md`` and the make files there for more information. Here we will give
a brief overview.

Besides building executable, other common make commands include:

    ``make clean``
        This removes the executable, .o files, and the temporarily generated
        files. Note that one can add additional targets to this rule using the
        double colon (::)

    ``make realclean``
        This removes all files generated by make.

    ``make help``
        This shows the rules for compilation.

    ``make print-xxx``
        This shows the value of variable xxx. This is very useful for debugging
        and tweaking the make system.

Compiler flags are set in ``amrex/Tools/GNUMake/comps/``. Note that variables
like ``CC`` and ``CFLAGS`` are reset in that directory and their values in
environment variables are disregarded.  Site-specific setups (e.g., the MPI
installation) are in ``amrex/Tools/GNUMake/sites/``, which includes a generic
setup in ``Make.unknown``. You can override the setup by having your own
``sites/Make.$(host_name)`` file, where variable ``host_name`` is your host
name in the make system and can be found via ``make print-host_name``.  You can
also have an ``amrex/Tools/GNUMake/Make.local`` file to override various
variables. See ``amrex/Tools/GNUMake/Make.local.template`` for more examples of
how to customize the build process.


.. _sec:build:local:

Specifying your own compiler
----------------------------

The ``amrex/Tools/GNUMake/Make.local`` file can also specify your own compile
commands by setting the valiables ``CXX``, ``CC``, ``FC``, and ``F90``. This
might be neccarry if your systems contains non-standard names for compiler
commands.

For example, the following ``amrex/Tools/GNUMake/Make.local`` builds AMReX
using a specific compiler (in this case ``gcc-8``) without MPI. Whenever
``USE_MPI``  is true, this confguration defaults to the appropriate
``mpixxx`` command:
::

    ifeq ($(USE_MPI),TRUE)
      CXX = mpicxx
      CC  = mpicc
      FC  = mpif90
      F90 = mpif90
    else
      CXX = g++-8
      CC  = gcc-8
      FC  = gfortran-8
      F90 = gfortran-8
    endif

For building with MPI, we assume ``mpicxx``, ``mpif90``, etc. provide access to
the correct underlying compilers.


.. _sec:build:macos:

GCC on macOS
------------

The example configuration above should also run on the latest macOS. On macOS
the default cxx compiler is clang, whereas the default fortran compiler is
gfortran. Sometimes it is good to avoid mixing compilers, in that case we can
use the ``Make.local`` to force using GCC. However, macOS' Xcode ships with its
own (woefully outdated) version of GCC (4.2.1). It is therefore recommended to
install GCC using the `homebrew <https://brew.sh>`_ package manager. Running
``brew install gcc`` installs gcc with names reflecting the version number. If
GCC 8.2 is installed, homebrew installs it as ``gcc-8``. AMReX can be built
using ``gcc-8`` (with and without MPI) by using the following
``amrex/Tools/GNUMake/Make.local``:

::

    CXX = g++-8
    CC  = gcc-8
    FC  = gfortran-8
    F90 = gfortran-8

    INCLUDE_LOCATIONS += /usr/local/include

The additional ``INCLUDE_LOCATIONS`` are installed using homebrew also. Note
that if you are building AMReX using homebrew's gcc, it is recommended that you
use homebrew's mpich. Normally is it fine to simply install its binaries:
``brew install mpich``. But if you are experiencing problems, we suggest
building mpich usinging homebrew's gcc: ``brew install mpich --cc=gcc-8``.


.. _sec:build:lib:

Building libamrex
=================

If an application code already has its own elaborated build system and wants to
use AMReX an external library, this might be your choice. In this approach, one
runs ``./configure``, followed by ``make`` and ``make install``.
Other make options include ``make distclean`` and ``make uninstall``.  In the top
AMReX directory, one can run ``./configure -h`` to show the various options for
the configure script. This approach is built on the AMReX GNU Make system. Thus
the section on :ref:`sec:build:make` is recommended if any fine tuning is
needed.  The result of ``./configure`` is ``GNUmakefile`` in the AMReX
top directory.  One can modify the make file for fine tuning.

.. _sec:build:cmake:

Building with CMake
===================

An alternative to the approach described in the section on :ref:`sec:build:lib`
is to install AMReX as an external library by using the CMake build system.  A
CMake build is a two-step process. First ``cmake`` is invoked to create
configuration files and makefiles in a chosen directory (``builddir``).  This
is roughly equivalent to running ``./configure`` (see the section on
:ref:`sec:build:lib`). Next, the actual build and installation are performed by
invoking ``make install`` from within builddir. This installs the library files
in a chosen installation directory (``installdir``).  If no installation path
is provided by the user, AMReX will be installed in /path/to/amrex/installdir.
The CMake build process is summarized as follows:

.. highlight:: console

::

    mkdir /path/to/builddir
    cd    /path/to/builddir
    cmake [options] -DCMAKE_INSTALL_PREFIX:PATH=/path/to/installdir  /path/to/amrex
    make  install

In the above snippet, ``[options]`` indicates one or more options for the
customization of the build, as described in the subsection on
:ref:`sec:build:cmake:options`.  Although the AMReX source could be used as
build directory, we advise against doing so.  After the installation is
complete, builddir can be removed.


Cmake and macOS
---------------

You can also specify your own compiler in cmake using the
``-DCMAKE_C_COMPILER`` and ``-DCMAKE_CXX_COMPILER`` options. While not strictly
necessary when using homebrew on macOS, it is highly recommended that the user
specifies ``-DCMAKE_C_COMPILER=$(which gcc-X) -DCMAKE_CXX_COMPILER=$(which
g++-X)`` (where X is the GCC version installed by homebrew) when using
gfortran. This is because homebrew's cmake defaults to the clang c/c++
compiler. Normaly clang plays well with gfortran, but if there are some issues,
we recommend telling cmake to use gcc for c/c++ also.


.. _sec:build:cmake:options:

Customization options
---------------------

AMReX configuration settings may be specified on the command line with the
``-D`` option.  For example, one can enable OpenMP support as follows:

.. highlight:: console

::

    cmake -DENABLE_OMP=1 -DCMAKE_INSTALL_PREFIX:PATH=/path/to/installdir  /path/to/amrex

The list of available option is reported in the table on :ref:`tab:cmakevar`
below.


.. raw:: latex

   \begin{center}

.. _tab:cmakevar:

.. table:: AMReX build options

   +------------------------------+-------------------------------------------------+-------------+-----------------+
   | Option Name                  | Description                                     | Default     | Possible values |
   +==============================+=================================================+=============+=================+
   | DEBUG                        |  Build AMReX in debug mode                      | OFF         | ONE, OFF        |
   +------------------------------+-------------------------------------------------+-------------+-----------------+
   | DIM                          |  Dimension of AMReX build                       | 3           | 2, 3            |
   +------------------------------+-------------------------------------------------+-------------+-----------------+
   | ENABLE_DP                    |  Build with double-precision reals              | ON          | ON, OFF         |
   +------------------------------+-------------------------------------------------+-------------+-----------------+
   | ENABLE_PIC                   |  Build Position Independent Code                | OFF         | ON, OFF         |
   +------------------------------+-------------------------------------------------+-------------+-----------------+
   | ENABLE_MPI                   |  Build with MPI support                         | ON          | ON OFF          |
   +------------------------------+-------------------------------------------------+-------------+-----------------+
   | ENABLE_OMP                   |  Build with OpenMP support                      | OFF         | ON, OFF         |
   +------------------------------+-------------------------------------------------+-------------+-----------------+
   | ENABLE_FORTRAN_INTERFACES    |  Build Fortran API                              | ON          | ON, OFF         |
   +------------------------------+-------------------------------------------------+-------------+-----------------+
   | ENABLE_LINEAR_SOLVERS        |  Build AMReX linear solvers                     | ON          | ON, OFF         |
   +------------------------------+-------------------------------------------------+-------------+-----------------+
   | ENABLE_LINEAR_SOLVERS_LEGACY |  Build AMReX linear solvers (legacy components) | ON          | ON, OFF         |
   +------------------------------+-------------------------------------------------+-------------+-----------------+
   | ENABLE_FBASELIB              |  Build (deprecated) Fortran kernel              | ON          | ON, OFF         |
   +------------------------------+-------------------------------------------------+-------------+-----------------+
   | ENABLE_AMRDATA               |  Build data services                            | OFF         | ON, OFF         |
   +------------------------------+-------------------------------------------------+-------------+-----------------+
   | ENABLE_PARTICLES             |  Build particle classes                         | OFF         | ON OFF          |
   +------------------------------+-------------------------------------------------+-------------+-----------------+
   | ENABLE_DP_PARTICLES          |  Use double-precision reals in particle classes | ON          | ON, OFF         |
   +------------------------------+-------------------------------------------------+-------------+-----------------+
   | ENABLE_BASE_PROFILE          |  Build with basic profiling support             | OFF         | ON, OFF         |
   +------------------------------+-------------------------------------------------+-------------+-----------------+
   | ENABLE_TINY_PROFILE          |  Build with tiny profiling support              | OFF         | ON, OFF         |
   +------------------------------+-------------------------------------------------+-------------+-----------------+
   | ENABLE_TRACE_PROFILE         |  Build with trace-profiling support             | OFF         | ON, OFF         |
   +------------------------------+-------------------------------------------------+-------------+-----------------+
   | ENABLE_COMM_PROFILE          |  Build with comm-profiling support              | OFF         | ON, OFF         |
   +------------------------------+-------------------------------------------------+-------------+-----------------+
   | ENABLE_MEM_PROFILE           |  Build with memory-profiling support            | OFF         | ON, OFF         |
   +------------------------------+-------------------------------------------------+-------------+-----------------+
   | ENABLE_PROFPARSER            |  Build with profile parser support              | OFF         | ON, OFF         |
   +------------------------------+-------------------------------------------------+-------------+-----------------+
   | ENABLE_BACKTRACE             |  Build with backtrace support                   | OFF         | ON, OFF         |
   +------------------------------+-------------------------------------------------+-------------+-----------------+
   | ENABLE_FPE                   |  Build with Floating Point Exceptions checks    | OFF         | ON,OFF          |
   +------------------------------+-------------------------------------------------+-------------+-----------------+
   | ENABLE_ASSERTIONS            |  Build with assertions turned on                | OFF         | ON,OFF          |
   +------------------------------+-------------------------------------------------+-------------+-----------------+
   | CMAKE_Fortran_FLAGS          |  User-defined Fortran flags                     |             | user-defined    |
   +------------------------------+-------------------------------------------------+-------------+-----------------+
   | CMAKE_CXX_FLAGS              |  User-defined C++ flags                         |             | user-defined    |
   +------------------------------+-------------------------------------------------+-------------+-----------------+
   | ENABLE_3D_NODAL_MGML         |  Enable 3D nodal projection                     | OFF         | ON,OFF          |
   +------------------------------+-------------------------------------------------+-------------+-----------------+
   | ALGOIM_INSTALL_DIR           |  Path to Algoim installation directory          |             | user-defined    |
   +------------------------------+-------------------------------------------------+-------------+-----------------+
   | BLITZ_INSTALL_DIR            |  Path to Blitz installation directory           |             | user-defined    |
   +------------------------------+-------------------------------------------------+-------------+-----------------+
.. raw:: latex

   \end{center}

The option ``ENABLE_LINEAR_SOLVERS=ON`` triggers the inclusion of C++-based
linear solvers in the build. Fortran-based linear solvers can be included as
well by providing the option ``ENABLE_FBASELIB=ON`` in addition to
``ENABLE_LINEAR_SOLVERS=ON``.

The option ``DEBUG=ON`` implies ``ENABLE_ASSERTION=ON``. In order to turn off
assertions in debug mode, ``ENABLE_ASSERTION=OFF`` must be set explicitly while
invoking CMake.

The options ``CMAKE_Fortran_FLAGS`` and ``CMAKE_CXX_FLAGS`` allow the user to
set his own compilation flags for Fortran and C++ source files respectively.
If ``CMAKE_Fortran_FLAGS``/ ``CMAKE_CXX_FLAGS`` are not set by the user,
they will be initialized with the value of the environmental variables ``FFLAGS``/
``CXXFLAGS``. If ``FFLAGS``/ ``CXXFLAGS`` are not defined in the environment,
``CMAKE_Fortran_FLAGS``/ ``CMAKE_CXX_FLAGS`` will be set to the AMReX default values
defined in  ``/path/to/amrex/Tools/CMake/AMReX_Compilers.cmake``.

The option ``ENABLE_3D_NODAL_MGML`` enables AMReX 3D nodal projection. This option requires
two external libraries: Blitz and Algoim. The user can provide the location of
both libraries via ``BLITZ_INSTALL_DIR`` and ``ALGOIM_INSTALL_DIR``. However, if one or both of these
options is not provided, AMReX will download and build Blitz and/or Algoim automatically.
It should be noted that AMReX 2D nodal projection does not require the use of external libraries.


.. _sec:build:cmake:config:

Importing AMReX into your CMake project
--------------------------------------------------

In order to import the AMReX library into your CMake project, you need
to include the following line in the appropriate CMakeLists.txt file:

.. highlight:: cmake

::

    find_package (AMReX 18 [REQUIRED] [HINTS /path/to/installdir/] )


In the above snippet, ``18`` refer to the mininum AMReX version supporting
the import feature discussed here.
Linking AMReX to any target defined in your CMake project is done by including
the following line in the appropriate CMakeLists.txt file

.. highlight:: cmake

::

    target_link_libraries ( <your-target-name>  AMReX::amrex )

The above snippet will take care of properly linking ``<your-target-name>``
to AMReX and to all the required transitive dependencies.
