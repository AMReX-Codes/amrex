.. role:: cpp(code)
   :language: c++

.. role:: fortran(code)
   :language: fortran


Compiling AMReX with Sundials (Version 5.1.0)
---------------------------------------------

The following steps describe how to compile an AMReX application with
SUNDIALS_5 support.

In order to use Sundials 5.1.0:

#. AMReX suggests using the Github mirror:
   https://github.com/LLNL/sundials/tree/v5.1.0

   ::

      #!/bin/bash
      set -e
      git clone https://github.com/LLNL/sundials
      cd sundials
      git checkout v5.1.0
      mkdir builddir instdir
      INSTALL_PREFIX=$(pwd)/instdir
      cd builddir
      cmake     \
      -DCMAKE_INSTALL_PREFIX=${INSTALL_PREFIX}     \
      -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON     \
      -DCMAKE_C_COMPILER=$(which gcc)     \
      -DCMAKE_CXX_COMPILER=$(which g++)     \
      -DCMAKE_CUDA_HOST_COMPILER=$(which g++)    \
      -DEXAMPLES_INSTALL_PATH=${INSTALL_PREFIX}/examples     \
      -DCMAKE_BUILD_TYPE=Release     \
      -DCMAKE_C_FLAGS_RELEASE="-O3 -DNDEBUG"     \
      -DCMAKE_CXX_FLAGS_RELEASE="-O3 -DNDEBUG"     \
      -DCUDA_ENABLE=ON     \
      -DMPI_ENABLE=OFF     \
      -DOPENMP_ENABLE=ON     \
      -DEXAMPLES_ENABLE=ON      \
      -DCUDA_ARCH=sm_70 ../
      make -j8
      make install -j8

#. Note that CMAKE_C_COMPILER and CMAKE_CXX_COMPILER need to be
   consistent with the AMReX make variable COMP to ensure
   matching OMP runtime libraries for use with the OpenMP NVector. 

#. For more detailed instructions for installing Sundials with different flags and versions see :ref:`CVODE`.

#. In the ``GNUmakefile`` for the  application which uses the Fortran 2003
   interface to CVODE or ARKODE, add ``USE_SUNDIALS = TRUE``, which will compile the Fortran 2003
   interfaces and link the  libraries.  Note that one must define the
   ``CVODE_LIB_DIR`` environment variable to point to the location where the
   libraries are installed. This should be an absolute path, such as ``$(pwd)/../sundials/instdir/lib64``.

#. In the ``GNUmakefile`` for the  application which uses the Fortran 2003
   interface to ARKODE, also add ``USE_ARKODE_LIBS = TRUE``. It is assumed that the
   ``CVODE_LIB_DIR`` environment variable points to the location where the ARKODE
   libraries are installed as well.

#. Fortran 2003 interfaces for the pgi compilers are currently not supported.

SUNDIALS 5 Tutorials
--------------------------

AMReX provides six tutorials in the ``amrex/Tutorials/CVODE/SUNDIALS3_finterface`` directory and
five tutorials in the ``amrex/Tutorials/SUNDIALS`` directory. See the Tutorials SUNDIALS_ documentation for more detail.

.. _SUNDIALS: https://amrex-codes.github.io/amrex/tutorials_html/SUNDIALS_Tutorial.html
