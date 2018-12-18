.. role:: cpp(code)
   :language: c++

.. role:: fortran(code)
   :language: fortran


Compiling AMReX with Sundials version 3.X or later
---------------------------------------------------

The following steps describe how to compile an AMReX application with
SUNDIALS_3.X support.  On Cray systems (e.g., Cori or Edison at NERSC), Cray provides
a system module called ``cray-tpsl`` (“Cray Third-Party Scientific Libraries”)
which as of this writing contains the 2.7 version of the SUNDIALS solver suite (including
CVODE).  

In order to use the Sundials 3.X version:

#. Obtain the CVODE source code, which is hosted here:
   https://computation.llnl.gov/projects/sundials/sundials-software.
   One can download either the complete SUNDIALS package, or just the CVODE components.

#. Unpack the CVODE / SUNDIALS tarball, and create a new “build” directory (it
   can be anywhere).

#. Navigate to the new, empty build directory, and type

   ::

         cmake \
           -DCMAKE_INSTALL_PREFIX:PATH=/path/to/install/dir \
           /path/to/cvode/or/sundials/top/level/source/dir


   The ``CMAKE_INSTALL_DIR`` option tells CMake where to install the libraries.
   Note that CMake will attempt to deduce the compilers automatically, but
   respects certain environment variables if they are defined, such as ``CC``
   (for the C compiler), ``CXX`` (for the C++ compiler), and ``FC`` (for the
   Fortran compiler).  So one may modify the above CMake invocation to be
   something like the following:

   ::

         CC=/path/to/gcc \
         CXX=/path/to/g++ \
         FC=/path/to/gfortran \
           cmake \
           -DCMAKE_INSTALL_PREFIX:PATH=/path/to/install/dir \
           /path/to/cvode/or/sundials/top/level/source/dir


   One can supply additional flags to CMake or to the compiler to customize the
   compilation process.  Flags of interest may include ``CMAKE_C_FLAGS``, which
   add the specified flags to the compile statement, e.g.,
   ``-DCMAKE_C_FLAGS="-h list=a"`` will append the ``-h list=a`` flag to the
   ``cc`` statement when compiling the source code.  Here one may wish to add
   something like ``"-O2 -g"`` to provide an optimized library that still
   contains debugging symbols; if one neglects debugging symbols in the CVODE
   library, and if a code that uses CVODE encounters a segmentation fault in
   the solve, then the backtrace has no information about where in the solver
   the error occurred.  Also, if one wishes to compile only the solver library
   itself and not the examples that come with the source (compiling the
   examples is enabled by default), one can add ``"-DEXAMPLES_ENABLE=OFF"``.
   Users should be aware that the CVODE examples are linked dynamically, so
   when compiling the solver library on Cray system using the Cray compiler
   wrappers ``cc``, ``CC``, and ``ftn``, one should explicitly disable
   compiling the examples via the ``"-DEXAMPLES_ENABLE=OFF"`` flag.

#. In the ``GNUmakefile`` for the  application which uses the Fortran 2003
   interface to CVODE or ARKODE, add ``SUNDIALS_3x4x = TRUE``, which will compile the Fortran 2003
   interfaces and link the  libraries.  Note that one must define the
   ``CVODE_LIB_DIR`` environment variable to point to the location where the
   libraries are installed.

#. In the ``GNUmakefile`` for the  application which uses the Fortran 2003
   interface to ARKODE, also add ``USE_ARKODE_LIBS = TRUE``. It is assumed that the
   ``CVODE_LIB_DIR`` environment variable points to the location where the ARKODE
   libraries are installed as well.

#. Fortran 2003 interfaces for the pgi compilers and for developmental versions of SUNDIALS
   are currently not supported.

SUNDIALS 3.X Tutorials
-------------------------

AMReX provides six tutorials in the ``amrex/Tutorials/CVODE/SUNDIALS3_finterface`` directory.
``EX1`` is modeled after the CVODE Tutorial ``EX1`` showing use with AMReX.
The four ``EX_cv_*`` tutorials are based on examples provided with the interface, which
are more closely modeled after CVODE examples. The ``EX_ark_analytic_fp`` tutorial is based
on the ``EX_cv_analytic_fp`` tutorial, but uses ARKODE instead of CVODE.

AMReX provides three tutorials in the ``amrex/Tutorials/CVODE/SUNDIALS3_cppversion`` directory.
These are versions of ``EX1`` which operate on a packed version of the data. ``EX1_SERIAL_NVEC``
packs a box worth of equations into a serial NVector, uses CVODE to solve, and then unpacks
the solution back into the box it came from. ``EX1_CUDA_NVEC`` uses the cuda NVector implementation instead.
``EX1_GPU_PRAGMA`` uses the cuda NVector, and the gpu pragma functionality.

.. _SUNDIALS3:
