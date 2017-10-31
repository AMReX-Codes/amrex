Introduction
============================

In this chapter, we discuss AMReX's build systems. There are three
ways to use AMReX. The approach used by AMReX developers uses GNU
Make. There is no installation step in this approach. Application
codes adopt AMReX's build system and compile AMReX while compiling
their own codes. This will be discussed in more detail in
Section 3.1. The second approach is to build into a library and 
install it (Section 3.2). Then an application code uses its own build system 
and links AMReX as an external library. AMReX can also be built with Cmake
(Section 3.3).

Building with GNU Make
======================

In this build approach, you write your own make files defining a
number of variables and rules. Then you invoke make to start
the building process. This will result in an executable upon
successful completion. The temporary files generated in the building
process are stored in a temporary directory named tmp\_build\_dir.

Dissecting a Simple Make File
-----------------------------

An example of building with GNU Make can be found in amrex/Tutorials/Basic/HelloWorld\_C. 
The table below shows a list of important variables.

+---------------+---------------------------------------+---------------+
| Variable      | Value                                 | Default       |
+===============+=======================================+===============+
| AMREX\_HOME   | Path to amrex                         | environment   |
+---------------+---------------------------------------+---------------+
| COMP          | gnu, cray, ibm, intel, llvm, or pgi   | none          |
+---------------+---------------------------------------+---------------+
| DEBUG         | TRUE or FALSE                         | TRUE          |
+---------------+---------------------------------------+---------------+
| DIM           | 1 or 2 or 3                           | none          |
+---------------+---------------------------------------+---------------+
| USE\_MPI      | TRUE or FALSE                         | FALSE         |
+---------------+---------------------------------------+---------------+
| USE\_OMP      | TRUE or FALSE                         | FALSE         |
+---------------+---------------------------------------+---------------+

Table: Important make variables

At the beginning of amrex/Tutorials/Basic/HelloWorld\_C/GNUmakefile, AMREX\_HOME
is set to the path to the top directory of AMReX. Note that in the
example ?= is a conditional variable assignment operator that
only has an effect if AMREX\_HOME has not been defined
(including in the environment). One can also set AMREX\_HOME
as an environment variable. For example in bash,
one can set export AMREX\_HOME=/path/to/amrex; in tcsh one can set
setenv AMREX\_HOME /path/to/amrex.

One must set the COMP variable to choose a compiler. Currently
the list of supported compilers includes gnu, cray, ibm, intel, llvm, and pgi. 
One must also set the DIM variable to either 1, 2, or 3, depending on the
dimensionality of the problem.

Variables DEBUG, USE\_MPI and USE\_OMP are optional
with default set to TRUE, FALSE and FALSE,
respectively. The meaning of these variables should be obvious.
When DEBUG = TRUE, aggressive compiler optimization flags are turned
off and assertions in AMReX source code are turned on. For
production runs, DEBUG should be set to FALSE.

After defining these make variables, a number of files, 
Make.defs, Make.package and Make.rules, are included in
GNUmakefile. based applications do not need to include
all directories in AMReX; an application which does not use particles,
for example, does not need to include files from the Particle
directory in its build.
In this simple example, we only need to include $(AMREX\_HOME)/Src/Base/Make.package. 
An application code also has its own Make.package file (e.g., ./Make.package in
this example) to append source files to the build system using
operator +=. Variables for various source files are shown
below.

    CEXE\_sources
         source files. Note that source files are assumed to have a .cpp extension.

    CEXE\_headers
         headers with .h or .H extension.

    cEXE\_sources
        C source files with .c extension.

    cEXE\_headers
        C headers with .h extension.

    f90EXE\_sources
        Free format Fortran source with .f90 extension.

    F90EXE\_sources
        Free format Fortran source with .F90 extension. Note that these Fortran files will go through
        preprocessing.

In this simple example, the extra source file, main.cpp is in
the current directory that is already in the build system's search
path. If this example has files in a subdirectory (e.g., mysrcdir), 
you will then need to add the following to Make.package.

::

        VPATH_LOCATIONS += mysrcdir
        INCLUDE_LOCATIONS += mysrcdir

Here VPATH\_LOCATIONS and INCLUDE\_LOCATIONS are the search
path for source and header files, respectively.

Tweaking Make System
--------------------

The GNU Make build system is located at amrex/Tools/GNUMake.
You can read README.md and the make files there for more
information. Here we will give a brief overview.

Besides building executable, other common make commands include:

    make clean
        This removes the executable, .o files, and
        the temporarily generated files. Note that one can add additional
        targets to this rule by using the double colon (::)

    make realclean
        This removes all files generated by make.

    make help
        This shows the rules for compilation.

    make print-xxx
        This shows the value of variable xxx. This is
        very useful for debugging and tweaking the make system.

Compiler flags are set in amrex/Tools/GNUMake/comps/. Note that
variables like CC and CFLAGS are reset in that directory
and their values in environment variables are disregarded. Site
specific setups (e.g., MPI installation) are in amrex/Tools/GNUMake/sites/, 
which includes a generic setup in Make.unknown. 
You can override the setup by having your own sites/Make.$(host\_name) file, 
where variable host\_name is your host name in the make system and can be found via make
print-host\_name. You can also have a amrex/Tools/GNUMake/Make.local file to override various variables.
See amrex/Tools/GNUMake/Make.local.template for an example.

Building libamrex
=================

If an application code already has its own elaborated build system and
wants to use AMReX as an external library, this might be your
choice. In this approach, one runs ./configure, followed by
make and make install. In the top AMReX directory, one
can run ./configure -h to show the various options for the configure script. 
This approach is built on the AMReX GNU Make
system. Thus Section 3.1 is recommended if any fine tuning is needed.

Building with CMake
===================

An alternative to the approach described in Section 3.2
is to install AMReX as an external library by using the CMake build system.
A CMake build is a two-steps process. First cmake is invoked to create
configuration files and makefiles in a chosen directory (builddir).
This is roughly equivalent to running ./configure (see Section 3.2).
Next, the actual build and installation are performed
by issuing make install from within builddir. This will install
the library files in a chosen installation directory ( installdir). 
If no installation path is provided by the user,
AMReX will be installed in /path/to/amrex/installdir.
The CMake build process is summarized as follow:

::

    mkdir /path/to/builddir
    cd    /path/to/builddir
    cmake [options] -DCMAKE_INSTALL_PREFIX:PATH=/path/to/installdir  /path/to/amrex 
    make install

In the above snippet, indicates one or more options for the customization
of the build, as described in Subsection 3.3.1.
Although the AMReX source could be used as build directory, we advise against doing so.
After the installation is complete, builddir can be removed.

Customization options
---------------------

 configuration settings may be specified on the command line with the -D option.
For example, one can enable OpenMP support as follows:

::

    cmake -DENABLE_OMP=1 -DCMAKE_INSTALL_PREFIX:PATH=/path/to/installdir  /path/to/amrex 

The list of available option is reported in the table below.

.. table:: Variables for the customization of AMReX build with CMake.

=========================== =============================================== ======== ===============
Option name                 Description                                     Defaults Possible values 
=========================== =============================================== ======== ===============
DEBUG                       Build \amrex\ in debug mode                     OFF      ON,OFF
DIM                         Dimension of \amrex\ build                      3        2,3
ENABLE\_DP                  Build with double-precision reals               ON       ON,OFF
ENABLE\_PIC                 Build Position Independent Code                 OFF      ON,OFF
ENABLE\_MPI                 Build with MPI support                          ON       ON,OFF
ENABLE\_OMP                 Build with OpenMP support                       OFF      ON,OFF
ENABLE\_FORTRAN\_INTERFACES Build Fortran API                               ON       ON,OFF
ENABLE\_LINEAR\_SOLVERS     Build \amrex\ linear solvers                    ON       ON,OFF
ENABLE\_FBASELIB            Build (deprecated) Fortran kernel               ON       ON,OFF
ENABLE\_AMRDATA             Build data services                             OFF      ON,OFF
ENABLE\_PARTICLES           Build particle classes                          OFF      ON,OFF
ENABLE\_DP\_PARTICLES       Use double-precision reals in particle classes  ON       ON,OFF
ENABLE\_BASE\_PROFILE       Build with basic profiling support              OFF      ON,OFF
ENABLE\_TINY\_PROFILE       Build with tiny profiling support               OFF      ON,OFF
ENABLE\_TRACE\_PROFILE      Build with trace-profiling support              OFF      ON,OFF
ENABLE\_COMM\_PROFILE       Build with comm-profiling support               OFF      ON,OFF
ENABLE\_MEM\_PROFILE        Build with memory-profiling support             OFF      ON,OFF
ENABLE\_PROFPARSER          Build with profile parser support               OFF      ON,OFF
ENABLE\_BACKTRACE           Build with backtrace support                    OFF      ON,OFF
ENABLE\_FPE                 Build with Floating Point Exceptions checks     OFF      ON,OFF
ENABLE\_ASSERTIONS          Build with assertions turned on                 OFF      ON,OFF
AMREX\_FFLAGS\_OVERRIDES    User-defined Fortran flags                      None     user-defined
AMREX\_CXXFLAGS\_OVERRIDES  User-defined C++ flags                          None     user-defined
=========================== =============================================== ======== ===============

The option ENABLE\_LINEAR\_SOLVERS=ON triggers the inclusion of C++-based linear
solvers in the build. Fortran-based linear solvers can be included
as well by providing the option ENABLE\_FBASELIB=ON in addition
to ENABLE\_LINEAR\_SOLVERS=ON.

The options DEBUG=ON implies ENABLE\_ASSERTION=ON. In order to turn off assertions in debug mode, 
ENABLE\_ASSERTION=OFF must be set explicitly while invoking CMake.

Importing AMReX configuration into a CMake project
--------------------------------------------------

In order to import the AMReX configuration options into your CMake
build system, include the following line in the appropriate
CMakeLists.txt file:

::

    find_package (AMReX CONFIG REQUIRED HINTS /path/to/installdir/cmake )

This will load -specific CMake variables containing the necessary
info to compile and link your code to . For a list of all the available
configuration variables, refer to the file AMReXConfig.cmake.in in
/path/to/installdir/cmake/.
