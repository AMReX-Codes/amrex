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

Note: when setting ``AMREX_HOME`` in the ``GNUmakefile``, be aware that ``~`` does
not expand, so ``AMREX_HOME=~/amrex/`` will yield an error.

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

If you need to pass macro definitions to the preprocessor, you can add
them to your make file as follows,

::

        DEFINES += -Dmyname1 -Dmyname2=mydefinition

To link to an additional library say ``foo`` with headers located at
``foopath/include`` and library at ``foopath/lib``, you can add the
following to your make file before the line that includes AMReX's
``Make.defs``,

::

        INCLUDE_LOCATIONS += foopath/include
        LIBRARY_LOCATIONS += foopath/lib
        LIBRARIES += -lfoo

.. _sec:build:local:

Specifying your own compiler
----------------------------

The ``amrex/Tools/GNUMake/Make.local`` file can also specify your own compile
commands by setting the variables ``CXX``, ``CC``, ``FC``, and ``F90``. This
might be necessary if your systems contains non-standard names for compiler
commands.

For example, the following ``amrex/Tools/GNUMake/Make.local`` builds AMReX
using a specific compiler (in this case ``gcc-8``) without MPI. Whenever
``USE_MPI``  is true, this configuration defaults to the appropriate
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
the default cxx compiler is clang, whereas the default Fortran compiler is
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
building mpich using homebrew's gcc: ``brew install mpich --cc=gcc-8``.


.. _sec:build:lib:

Building libamrex
=================

If an application code already has its own elaborated build system and wants to
use AMReX, an external AMReX library can be created instead. In this approach, one
runs ``./configure``, followed by ``make`` and ``make install``.
Other make options include ``make distclean`` and ``make uninstall``.  In the top
AMReX directory, one can run ``./configure -h`` to show the various options for
the configure script. In particular, one can specify the installation path for the AMReX library using::

  ./configure --prefix=[AMReX library path]

This approach is built on the AMReX GNU Make system. Thus
the section on :ref:`sec:build:make` is recommended if any fine tuning is
needed.  The result of ``./configure`` is ``GNUmakefile`` in the AMReX
top directory.  One can modify the make file for fine tuning.

To compile an application code against the external AMReX library, it
is necessary to set appropriate compiler flags and set the library
paths for linking. To assist with this, when the AMReX library is
built, a configuration file is created in ``[AMReX library path]/lib/pkgconfig/amrex.pc``.
This file contains the Fortran and
C++ flags used to compile the AMReX library as well as the appropriate
library and include entries.

The following sample GNU Makefile will compile a ``main.cpp`` source
file against an external AMReX library, using the C++ flags and
library paths used to build AMReX::

  AMREX_LIBRARY_HOME ?= [AMReX library path]

  LIBDIR := $(AMREX_LIBRARY_HOME)/lib
  INCDIR := $(AMREX_LIBRARY_HOME)/include

  COMPILE_CPP_FLAGS ?= $(shell awk '/Cflags:/ {$$1=$$2=""; print $$0}' $(LIBDIR)/pkgconfig/amrex.pc)
  COMPILE_LIB_FLAGS ?= $(shell awk '/Libs:/ {$$1=$$2=""; print $$0}' $(LIBDIR)/pkgconfig/amrex.pc)

  CFLAGS := -I$(INCDIR) $(COMPILE_CPP_FLAGS)
  LFLAGS := -L$(LIBDIR) $(COMPILE_LIB_FLAGS)

  all:
          g++ -o main.exe main.cpp $(CFLAGS) $(LFLAGS)

.. _sec:build:cmake:

Building with CMake
===================

An alternative to the approach described in the section on :ref:`sec:build:lib`
is to install AMReX as an external library by using the CMake build system.  A
CMake build is a two-step process. First ``cmake`` is invoked to create
configuration files and makefiles in a chosen directory (``builddir``).  This
is roughly equivalent to running ``./configure`` (see the section on
:ref:`sec:build:lib`). Next, the actual build and installation are performed by
invoking ``make install`` from within ``builddir``. This installs the library files
in a chosen installation directory (``installdir``).  If no installation path
is provided by the user, AMReX will be installed in ``/path/to/amrex/installdir``.
The CMake build process is summarized as follows:

.. highlight:: console

::

    mkdir /path/to/builddir
    cd    /path/to/builddir
    cmake [options] -DCMAKE_BUILD_TYPE=[Debug|Release|RelWithDebInfo|MinSizeRel] -DCMAKE_INSTALL_PREFIX=/path/to/installdir  /path/to/amrex
    make  install

In the above snippet, ``[options]`` indicates one or more options for the
customization of the build, as described in the subsection on
:ref:`sec:build:cmake:options`. If the option ``CMAKE_BUILD_TYPE`` is omitted,
``CMAKE_BUILD_TYPE=Release`` is assumed. Although the AMReX source could be used as
build directory, we advise against doing so.  After the installation is
complete, ``builddir`` can be removed.


.. _sec:build:cmake:options:

Customization options
---------------------

AMReX build can be customized  by setting the value of suitable configuration variables
on the command line via the ``-D <var>=<value>`` syntax, where ``<var>`` is the
variable to set and ``<value>`` its desired value.
For example, one can enable OpenMP support as follows:

.. highlight:: console

::

    cmake -DENABLE_OMP=YES -DCMAKE_INSTALL_PREFIX=/path/to/installdir  /path/to/amrex

In the example above ``<var>=ENABLE_OMP`` and ``<value>=YES``.
Configuration variables requiring a boolen value are evaluated to true if they
are assigned a value of ``1``, ``ON``, ``YES``, ``TRUE``, ``Y``. Conversely they are evaluated to false
if they are assigned a value of ``0``, ``OFF``, ``NO``, ``FALSE``, ``N``.
Boolean configuration variables are case-insensitive.
The list of available options is reported in the table on :ref:`tab:cmakevar`
below.


.. raw:: latex

   \begin{center}

.. _tab:cmakevar:

.. table:: AMReX build options

   +------------------------------+-------------------------------------------------+-------------+-----------------+
   | Variable Name                | Description                                     | Default     | Possible values |
   +==============================+=================================================+=============+=================+
   | CMAKE_Fortran_COMPILER       |  User-defined Fortran compiler                  |             | user-defined    |
   +------------------------------+-------------------------------------------------+-------------+-----------------+
   | CMAKE_CXX_COMPILER           |  User-defined C++ compiler                      |             | user-defined    |
   +------------------------------+-------------------------------------------------+-------------+-----------------+   
   | CMAKE_Fortran_FLAGS          |  User-defined Fortran flags                     |             | user-defined    |
   +------------------------------+-------------------------------------------------+-------------+-----------------+
   | CMAKE_CXX_FLAGS              |  User-defined C++ flags                         |             | user-defined    |
   +------------------------------+-------------------------------------------------+-------------+-----------------+
   | DIM                          |  Dimension of AMReX build                       | 3           | 1, 2, 3         |
   +------------------------------+-------------------------------------------------+-------------+-----------------+
   | USE_XSDK_DEFAULTS            |  Use XSDK defaults settings                     | NO          | YES, NO         |
   +------------------------------+-------------------------------------------------+-------------+-----------------+
   | ENABLE_DP                    |  Build with double-precision reals              | YES         | YES, NO         |
   +------------------------------+-------------------------------------------------+-------------+-----------------+
   | ENABLE_PIC                   |  Build Position Independent Code                | NO          | YES, NO         |
   +------------------------------+-------------------------------------------------+-------------+-----------------+
   | ENABLE_MPI                   |  Build with MPI support                         | YES         | YES, NO         |
   +------------------------------+-------------------------------------------------+-------------+-----------------+
   | ENABLE_OMP                   |  Build with OpenMP support                      | NO          | YES, NO         |
   +------------------------------+-------------------------------------------------+-------------+-----------------+
   | ENABLE_CUDA                  |  Build with CUDA support                        | NO          | YES, NO         |
   +------------------------------+-------------------------------------------------+-------------+-----------------+
   | CUDA_ARCH                    |  CUDA target architecture                       | Auto        | User-defined    |
   +------------------------------+-------------------------------------------------+-------------+-----------------+
   | CUDA_MAX_THREADS             |  Max number of CUDA threads per block           | 256         | User-defined    |
   +------------------------------+-------------------------------------------------+-------------+-----------------+
   | CUDA_MAXREGCOUNT             |  Limits the number of CUDA registers available  | 255         | User-defined    |
   +------------------------------+-------------------------------------------------+-------------+-----------------+
   | ENABLE_CUDA_FASTMATH         |  Enable CUDA fastmath library                   | YES         | YES, NO         |
   +------------------------------+-------------------------------------------------+-------------+-----------------+
   | ENABLE_FORTRAN_INTERFACES    |  Build Fortran API                              | NO          | YES, NO         |
   +------------------------------+-------------------------------------------------+-------------+-----------------+
   | ENABLE_LINEAR_SOLVERS        |  Build AMReX linear solvers                     | YES         | YES, NO         |
   +------------------------------+-------------------------------------------------+-------------+-----------------+
   | ENABLE_AMRDATA               |  Build data services                            | NO          | YES, NO         |
   +------------------------------+-------------------------------------------------+-------------+-----------------+
   | ENABLE_EB                    |  Build Embedded Boundary support                | NO          | YES, NO         |   
   +------------------------------+-------------------------------------------------+-------------+-----------------+
   | ENABLE_PARTICLES             |  Build particle classes                         | NO          | YES, NO         |
   +------------------------------+-------------------------------------------------+-------------+-----------------+
   | ENABLE_DP_PARTICLES          |  Use double-precision reals in particle classes | YES         | YES, NO         |
   +------------------------------+-------------------------------------------------+-------------+-----------------+
   | ENABLE_BASE_PROFILE          |  Build with basic profiling support             | NO          | YES, NO         |
   +------------------------------+-------------------------------------------------+-------------+-----------------+
   | ENABLE_TINY_PROFILE          |  Build with tiny profiling support              | NO          | YES, NO         |
   +------------------------------+-------------------------------------------------+-------------+-----------------+
   | ENABLE_TRACE_PROFILE         |  Build with trace-profiling support             | NO          | YES, NO         |
   +------------------------------+-------------------------------------------------+-------------+-----------------+
   | ENABLE_COMM_PROFILE          |  Build with comm-profiling support              | NO          | YES, NO         |
   +------------------------------+-------------------------------------------------+-------------+-----------------+
   | ENABLE_MEM_PROFILE           |  Build with memory-profiling support            | NO          | YES, NO         |
   +------------------------------+-------------------------------------------------+-------------+-----------------+
   | ENABLE_PROFPARSER            |  Build with profile parser support              | NO          | YES, NO         |
   +------------------------------+-------------------------------------------------+-------------+-----------------+
   | ENABLE_BACKTRACE             |  Build with backtrace support                   | NO          | YES, NO         |
   +------------------------------+-------------------------------------------------+-------------+-----------------+
   | ENABLE_FPE                   |  Build with Floating Point Exceptions checks    | NO          | YES, NO         |
   +------------------------------+-------------------------------------------------+-------------+-----------------+
   | ENABLE_ASSERTIONS            |  Build with assertions turned on                | NO          | YES, NO         |
   +------------------------------+-------------------------------------------------+-------------+-----------------+
   | ENABLE_SUNDIALS              |  Enable SUNDIALS 4 interfaces                   | NO          | YES, NO         |
   +------------------------------+-------------------------------------------------+-------------+-----------------+
   | ENABLE_SENSEI_IN_SITU        |  Enable SENSEI_IN_SITU infrastucture            | NO          | YES, NO         |
   +------------------------------+-------------------------------------------------+-------------+-----------------+
   | ENABLE_CONDUIT               |  Enable Conduit support                         | NO          | YES, NO         |
   +------------------------------+-------------------------------------------------+-------------+-----------------+
   | ENABLE_ASCENT                |  Enable Ascent support                          | NO          | YES, NO         |
   +------------------------------+-------------------------------------------------+-------------+-----------------+
   | ENABLE_HYPRE                 |  Enable HYPRE interfaces                        | NO          | YES, NO         |
   +------------------------------+-------------------------------------------------+-------------+-----------------+
   | ENABLE_PLOTFILE_TOOLS        |  Build and install plotfile postprocessing tools| NO          | YES, NO         |
   +------------------------------+-------------------------------------------------+-------------+-----------------+
   | ENABLE_TUTORIALS             |  Build tutorials                                | NO          | YES, NO         |
   +------------------------------+-------------------------------------------------+-------------+-----------------+   
.. raw:: latex

   \end{center}

The option ``CMAKE_BUILD_TYPE=Debug`` implies ``ENABLE_ASSERTION=YES``. In order to turn off
assertions in debug mode, ``ENABLE_ASSERTION=NO`` must be set explicitly while
invoking CMake.


The ``CMAKE_C_COMPILER``, ``CMAKE_CXX_COMPILER``, and  ``CMAKE_Fortran_COMPILER`` options
are used to tell CMake which compiler to use for the compilation of C, C++, and Fortran sources
respectively. If those options are not set by the user, CMake will use the system default compilers.

The options ``CMAKE_Fortran_FLAGS`` and ``CMAKE_CXX_FLAGS`` allow the user to
set his own compilation flags for Fortran and C++ source files respectively.
If ``CMAKE_Fortran_FLAGS``/ ``CMAKE_CXX_FLAGS`` are not set by the user,
they will be initialized with the value of the environmental variables ``FFLAGS``/
``CXXFLAGS``. If neither ``FFLAGS``/ ``CXXFLAGS`` nor ``CMAKE_Fortran_FLAGS``/ ``CMAKE_CXX_FLAGS``
are defined, AMReX default flags are used.

For a detailed explanation of CUDA support in AMReX CMake, refer to section :ref:`sec:gpu:build`.



CMake and macOS
---------------

While not strictly necessary when using homebrew on macOS, it is highly
recommended that the user specifies ``-DCMAKE_C_COMPILER=$(which gcc-X) -DCMAKE_CXX_COMPILER=$(which
g++-X)`` (where X is the GCC version installed by homebrew) when using
gfortran. This is because homebrew's CMake defaults to the Clang C/C++
compiler. Normally Clang plays well with gfortran, but if there are some issues,
we recommend telling CMake to use gcc for C/C++ also.

.. _sec:build:cmake:config:

Importing AMReX into your CMake project
--------------------------------------------------

In order to import AMReX into your CMake project, you need
to include the following line in the appropriate CMakeLists.txt file:

.. highlight:: cmake

::

    find_package(AMReX)


Calls to ``find_package(AMReX)`` will find a valid installation of AMReX, if present,
and import its settings and targets into your CMake project.
Imported AMReX targets can be linked to any of your targets, after they have been made available
following a successful call to ``find_package(AMReX)``, by including
the following line in the appropriate CMakeLists.txt file:

.. highlight:: cmake

::

    target_link_libraries( <your-target-name>  AMReX::<amrex-target-name> )

    
In the above snippet, ``<amrex-target-name>`` is any of the targets listed in the table below.

.. raw:: latex

   \begin{center}

.. _tab:cmaketargets:

.. table:: AMReX targets available for import. 

   +-----------------------+-------------------------------------------------+
   | Target name           | Description                                     |
   +=======================+=================================================+
   | amrex                 |  AMReX library                                  |
   +-----------------------+-------------------------------------------------+
   | Flags_CXX             |  C++ flags preset (interface)                   |
   +-----------------------+-------------------------------------------------+
   | Flags_Fortran         |  Fortran flags preset (interface)               |
   +-----------------------+-------------------------------------------------+
   | Flags_FPE             |  Floating Point Exception flags (interface)     |
   +-----------------------+-------------------------------------------------+
.. raw:: latex

   \end{center}


The options used to configure the AMReX build may result in certain parts, or ``components``, of the AMReX source code
to be excluded from compilation. For example, setting ``-DENABLE_LINEAR_SOLVERS=no`` at configure time
prevents the compilation of AMReX linear solvers code. 
Your CMake project can check which component is included in the AMReX library via `find_package`:


.. highlight:: cmake

::

    find_package(AMReX REQUIRED <components-list>)


The keyword ``REQUIRED`` in the snippet above will cause a fatal error if AMReX is not found, or
if it is found but the components listed in ``<components-list>`` are not include in the installation.
A list of AMReX component names and related configure options are shown in the table below.
  

.. raw:: latex

   \begin{center}

.. _tab:cmakecomponents:

.. table:: AMReX components.

   +------------------------------+-----------------+
   | Option                       | Component       |
   +==============================+=================+
   | DIM                          | 1D, 2D, 3D      |
   +------------------------------+-----------------+
   | ENABLE_DP                    | DP              |
   +------------------------------+-----------------+
   | ENABLE_PIC                   | PIC             |
   +------------------------------+-----------------+
   | ENABLE_MPI                   | MPI             |
   +------------------------------+-----------------+
   | ENABLE_OMP                   | OMP             |
   +------------------------------+-----------------+
   | ENABLE_CUDA                  | CUDA            |
   +------------------------------+-----------------+
   | ENABLE_FORTRAN_INTERFACES    | FINTERFACES     |
   +------------------------------+-----------------+
   | ENABLE_LINEAR_SOLVERS        | LSOLVERS        |
   +------------------------------+-----------------+
   | ENABLE_AMRDATA               | AMRDATA         |
   +------------------------------+-----------------+
   | ENABLE_EB                    | EB              |
   +------------------------------+-----------------+
   | ENABLE_PARTICLES             | PARTICLES       |
   +------------------------------+-----------------+
   | ENABLE_DP_PARTICLES          | DPARTICLES      |
   +------------------------------+-----------------+
   | ENABLE_BASE_PROFILE          | BASEP           |
   +------------------------------+-----------------+
   | ENABLE_TINY_PROFILE          | TINYP           |
   +------------------------------+-----------------+
   | ENABLE_TRACE_PROFILE         | TRACEP          |
   +------------------------------+-----------------+
   | ENABLE_COMM_PROFILE          | COMMP           |
   +------------------------------+-----------------+
   | ENABLE_MEM_PROFILE           | MEMP            |
   +------------------------------+-----------------+
   | ENABLE_PROFPARSER            | PROFPARSER      |
   +------------------------------+-----------------+
   | ENABLE_BACKTRACE             | BACKTRACE       |
   +------------------------------+-----------------+
   | ENABLE_FPE                   | FPE             |
   +------------------------------+-----------------+
   | ENABLE_ASSERTIONS            | ASSERTIONS      |
   +------------------------------+-----------------+
   | ENABLE_SUNDIALS              | SUNDIALS        |
   +------------------------------+-----------------+
   | ENABLE_SENSEI_IN_SITU        | SENSEI          |
   +------------------------------+-----------------+
   | ENABLE_CONDUIT               | CONDUIT         |
   +------------------------------+-----------------+
   | ENABLE_ASCENT                | ASCENT          |
   +------------------------------+-----------------+
   | ENABLE_HYPRE                 | HYPRE           |
   +------------------------------+-----------------+

.. raw:: latex

   \end{center}
   
As an example, consider the following CMake code:


.. highlight:: cmake

::

    find_package(AMReX REQUIRED 3D EB)
    target_link_libraries( Foo  AMReX::amrex AMReX::Flags_CXX )

The code in the snippet above checks wheather an AMReX installation with 3D and Embedded Boundary support
is available on the system. If so, AMReX is linked to target ``Foo`` and AMReX flags preset is used
to compile ``Foo``'s C++ sources. If no AMReX installation is found or if the available one was built without
3D or Embedded Boundary support, a fatal error is issued.


..
   It will fail if
   it cannot find any, or if the available one was not built with 3D and Embedded Boudary support.
   If AMReX is found, it will then link AMReX to target ``Foo`` and use the AMReX flags preset
   to compile ``Foo``'s C++ sources.


You can tell CMake to look for the AMReX library in non-standard paths by setting the environmental variable
``AMReX_ROOT`` to point to the AMReX installation directory or by adding
``-DAMReX_ROOT=<path/to/amrex/installation/directory>`` to the ``cmake`` invocation.
More details on ``find_package`` can be found 
`here <https://cmake.org/cmake/help/v3.14/command/find_package.html>`_.

