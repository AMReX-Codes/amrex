.. _sec:build:make:

Building with GNU Make
======================

In this build approach, you write your own make files defining a
number of variables and rules. Then you invoke make to start
the building process. This will result in an executable upon
successful completion. The temporary files generated in the building
process are stored in a temporary directory named tmp_build_dir.

Dissecting a Simple Make File
-----------------------------

An example of building with GNU Make can be found in
amrex/Tutorials/Basic/HelloWorld_C. The table on :ref:`tab:makevarimp`
below shows a list of important variables.

.. raw:: latex

   \centering

.. _tab:makevarimp:
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

At the beginning of amrex/Tutorials/Basic/HelloWorld_C/GNUmakefile, 
``AMREX_HOME`` is set to the path to the top directory of AMReX. 
Note that in the example ``?=`` is a conditional variable assignment 
operator that only has an effect if ``AMREX_HOME`` has not been defined
(including in the environment). One can also set ``AMREX_HOME`` as an 
environment variable. For example in bash, one can set 

::

    export AMREX_HOME=/path/to/amrex

alternatively, in tcsh one can set

::

    setenv AMREX_HOME /path/to/amrex

One must set the ``COMP`` variable to choose a compiler. Currently
the list of supported compilers includes gnu, cray,
ibm, intel, llvm, and pgi. One must also set the
DIM variable to either 1, 2, or 3, depending on the
dimensionality of the problem.

Variables ``DEBUG``, ``USE_MPI`` and ``USE_OMP`` are optional
with default set to TRUE, FALSE and FALSE, respectively. 
The meaning of these variables should be obvious.
When ``DEBUG = TRUE``, aggressive compiler optimization flags are turned
off and assertions in  source code are turned on. For
production runs, ``DEBUG`` should be set to FALSE.

After defining these make variables, a number of files,
Make.defs, Make.package and Make.rules, are included in
GNUmakefile. AMReX-based applications do not need to include
all directories in AMReX; an application which does not use particles,
for example, does not need to include files from the Particle
directory in its build.
In this simple example, we only need to include
``$(AMREX_HOME)/Src/Base/Make.package``. An application code also
has its own Make.package file (e.g., ./Make.package in
this example) to append source files to the build system using
operator ``+=``. Variables for various source files are shown
below.

    CEXE_sources
        C++ source files. Note that C++ source files are assumed to have a .cpp extension.

    CEXE_headers
        C++ headers with .h or .H extension.

    cEXE_sources
        C source files with .c extension.

    cEXE_headers
        C headers with .h extension.

    f90EXE_sources
        Free format Fortran source with .f90 extension.

    F90EXE_sources
        Free format Fortran source with .F90 extension. 
        Note that these Fortran files will go through preprocessing.

In this simple example, the extra source file, main.cpp is in
the current directory that is already in the build system's search
path. If this example has files in a subdirectory (e.g.,
mysrcdir), you will then need to add the following to
Make.package.

::

        VPATH_LOCATIONS += mysrcdir
        INCLUDE_LOCATIONS += mysrcdir

Here ``VPATH_LOCATIONS`` and ``INCLUDE_LOCATIONS`` are the search
path for source and header files, respectively.

Tweaking Make System
--------------------

The GNU Make build system is located at amrex/Tools/GNUMake.
You can read README.md and the make files there for more
information. Here we will give a brief overview.

Besides building executable, other common make commands include:

    ``make clean``
        This removes the executable, .o files, and
        the temporarily generated files. Note that one can add additional
        targets to this rule by using the double colon (::)

    ``make realclean``
        This removes all files generated by make.

    ``make help``
        This shows the rules for compilation.

    ``make print-xxx``
        This shows the value of variable xxx. This is
        very useful for debugging and tweaking the make system.

Compiler flags are set in amrex/Tools/GNUMake/comps/. Note that
variables like ``CC`` and ``CFLAGS`` are reset in that directory
and their values in environment variables are disregarded. Site
specific setups (e.g., MPI installation) are in
amrex/Tools/GNUMake/sites/, which includes a generic setup in
Make.unknown. You can override the setup by having your own
``sites/Make.$(host_name)`` file, where variable ``host_name`` is your
host name in the make system and can be found via ``make print-host_name``. 
You can also have an amrex/Tools/GNUMake/Make.local file to override 
various variables. See amrex/Tools/GNUMake/Make.local.template for an example.

.. _sec:build:lib:

Building libamrex
=================

If an application code already has its own elaborated build system and
wants to use  as an external library, this might be your
choice. In this approach, one runs ``./configure``, followed by
``make`` and ``make install``. In the top AMReX directory, one
can run ``./configure -h`` to show the various options for the
configure script. This approach is built on the AMReX GNU Make
system. Thus the section on :ref:`sec:build:make` is recommended if any fine
tuning is needed.

.. _sec:build:cmake:

Building with CMake
===================

An alternative to the approach described in the section on :ref:`sec:build:lib`
is to install AMReX as an external library by using the CMake build system.
A CMake build is a two-steps process. First ``cmake`` is invoked to create
configuration files and makefiles in a chosen directory (builddir).
This is roughly equivalent to running ``./configure`` (see the section on
:ref:`sec:build:lib`). Next, the actual build and installation are performed
by issuing ``make install`` from within builddir. This will install
the library files in a chosen installation directory (installdir). 
If no installation path is provided by the user, AMReX will be installed in 
/path/to/amrex/installdir. The CMake build process is summarized as follow:

::

    mkdir /path/to/builddir
    cd    /path/to/builddir
    cmake [options] -DCMAKE_INSTALL_PREFIX:PATH=/path/to/installdir  /path/to/amrex 
    make  install

In the above snippet, ``[options]`` indicates one or more options for the customization
of the build, as described in the subsection on :ref:`sec:build:cmake:options`.
Although the AMReX source could be used as build directory, we advise against doing so.
After the installation is complete, builddir can be removed.

.. _sec:build:cmake:options:

Customization options
---------------------

AMReX configuration settings may be specified on the command line with the -D option.
For example, one can enable OpenMP support as follows:

::

    cmake -DENABLE_OMP=1 -DCMAKE_INSTALL_PREFIX:PATH=/path/to/installdir  /path/to/amrex 

The list of available option is reported in the table on :ref:`tab:cmakevar` below.


.. raw:: latex

   \centering

.. _tab:cmakevar:
.. table:: Important cmake build options

   +---------------------------+-------------------------------------------------+-------------+-----------------+
   | Option Name               | Description                                     | Default     | Possible values |
   +===========================+=================================================+=============+=================+
   | DEBUG                     |  Build AMReX in debug mode                      | OFF         | ONE, OFF        |
   +---------------------------+-------------------------------------------------+-------------+-----------------+
   | DIM                       |  Dimension of AMReX build                       | 3           | 2, 3            |
   +---------------------------+-------------------------------------------------+-------------+-----------------+
   | ENABLE_DP                 |  Build with double-precision reals              | ON          | ON, OFF         |
   +---------------------------+-------------------------------------------------+-------------+-----------------+
   | ENABLE_PIC                |  Build Position Independent Code                | OFF         | ON, OFF         |
   +---------------------------+-------------------------------------------------+-------------+-----------------+
   | ENABLE_MPI                |  Build with MPI support                         | ON          | ON OFF          |
   +---------------------------+-------------------------------------------------+-------------+-----------------+
   | ENABLE_OMP                |  Build with OpenMP support                      | OFF         | ON, OFF         |
   +---------------------------+-------------------------------------------------+-------------+-----------------+
   | ENABLE_FORTRAN_INTERFACES |  Build Fortran API                              | ON          | ON, OFF         |
   +---------------------------+-------------------------------------------------+-------------+-----------------+
   | ENABLE_LINEAR_SOLVERS     |  Build AMReX linear solvers                     | ON          | ON, OFF         |
   +---------------------------+-------------------------------------------------+-------------+-----------------+
   | ENABLE_FBASELIB           |  Build (deprecated) Fortran kernel              | ON          | ON, OFF         |
   +---------------------------+-------------------------------------------------+-------------+-----------------+
   | ENABLE_AMRDATA            |  Build data services                            | OFF         | ON, OFF         |
   +---------------------------+-------------------------------------------------+-------------+-----------------+
   | ENABLE_PARTICLES          |  Build particle classes                         | OFF         | ON OFF          |
   +---------------------------+-------------------------------------------------+-------------+-----------------+
   | ENABLE_DP_PARTICLES       |  Use double-precision reals in particle classes | ON          | ON, OFF         |
   +---------------------------+-------------------------------------------------+-------------+-----------------+
   | ENABLE_BASE_PROFILE       |  Build with basic profiling support             | OFF         | ON, OFF         |
   +---------------------------+-------------------------------------------------+-------------+-----------------+
   | ENABLE_TINY_PROFILE       |  Build with tiny profiling support              | OFF         | ON, OFF         |
   +---------------------------+-------------------------------------------------+-------------+-----------------+
   | ENABLE_TRACE_PROFILE      |  Build with trace-profiling support             | OFF         | ON, OFF         |
   +---------------------------+-------------------------------------------------+-------------+-----------------+
   | ENABLE_COMM_PROFILE       |  Build with comm-profiling support              | OFF         | ON, OFF         |
   +---------------------------+-------------------------------------------------+-------------+-----------------+
   | ENABLE_MEM_PROFILE        |  Build with memory-profiling support            | OFF         | ON, OFF         | 
   +---------------------------+-------------------------------------------------+-------------+-----------------+
   | ENABLE_PROFPARSER         |  Build with profile parser support              | OFF         | ON, OFF         |
   +---------------------------+-------------------------------------------------+-------------+-----------------+
   | ENABLE_BACKTRACE          |  Build with backtrace support                   | OFF         | ON, OFF         |
   +---------------------------+-------------------------------------------------+-------------+-----------------+
   | ENABLE_FPE                |  Build with Floating Point Exceptions checks    | OFF         | ON,OFF          |
   +---------------------------+-------------------------------------------------+-------------+-----------------+
   | ENABLE_ASSERTIONS         |  Build with assertions turned on                | OFF         | ON,OFF          |
   +---------------------------+-------------------------------------------------+-------------+-----------------+
   | AMREX_FFLAGS_OVERRIDES}   |  User-defined Fortran flags                     | None        | user-defined    |
   +---------------------------+-------------------------------------------------+-------------+-----------------+
   | AMREX_CXXFLAGS_OVERRIDES} |  User-defined C++ flags                         | None        | user-defined    |
   +---------------------------+-------------------------------------------------+-------------+-----------------+


The option ``ENABLE_LINEAR_SOLVERS=ON`` triggers the inclusion of C++-based linear
solvers in the build. Fortran-based linear solvers can be included as well by providing 
the option ``ENABLE_FBASELIB=ON`` in addition to ``ENABLE_LINEAR_SOLVERS=ON``. The 
options ``DEBUG=ON`` implies ``ENABLE_ASSERTION=ON``. In order to turn off assertions 
in debug mode, ``ENABLE_ASSERTION=OFF`` must be set explicitly while invoking CMake.

.. _sec:build:cmake:config:

Importing AMReX configuration into a CMake project
--------------------------------------------------

In order to import the AMReX configuration options into your CMake
build system, include the following line in the appropriate
CMakeLists.txt file:

::

    find_package (AMReX CONFIG REQUIRED HINTS /path/to/installdir/cmake )

This will load AMReX-specific CMake variables containing the necessary
info to compile and link your code to AMReX. For a list of all the available
configuration variables, refer to the file AMReXConfig.cmake.in in
/path/to/installdir/cmake/.
