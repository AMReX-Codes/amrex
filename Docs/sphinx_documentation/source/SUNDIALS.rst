.. role:: cpp(code)
   :language: c++

.. role:: fortran(code)
   :language: fortran


Compiling AMReX with SUNDIALS 5
===============================

The following steps describe how to compile an AMReX application with
SUNDIALS 5 support.

In order to use SUNDIALS:

#. AMReX suggests using the Github mirror:
   https://github.com/LLNL/sundials

   ::

      #!/bin/bash
      set -e
      git clone https://github.com/LLNL/sundials
      cd sundials
      mkdir builddir instdir
      INSTALL_PREFIX=$(pwd)/instdir
      cd builddir
      cmake \
      -DCMAKE_INSTALL_PREFIX=${INSTALL_PREFIX}  \
      -DCMAKE_INSTALL_LIBDIR=lib \
      -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON \
      -DCMAKE_C_COMPILER=$(which gcc)  \
      -DCMAKE_CXX_COMPILER=$(which g++)   \
      -DCMAKE_CUDA_HOST_COMPILER=$(which g++)    \
      -DEXAMPLES_INSTALL_PATH=${INSTALL_PREFIX}/examples \
      -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_C_FLAGS_RELEASE="-O3 -DNDEBUG" \
      -DCMAKE_CXX_FLAGS_RELEASE="-O3 -DNDEBUG"  \
      -DCUDA_ENABLE=ON  \
      -DMPI_ENABLE=OFF  \
      -DOPENMP_ENABLE=ON   \
      -DF2003_INTERFACE_ENABLE=ON   \
      -DCUDA_ARCH=sm_70 ../
      make -j8
      make install -j8

#. Note that ``CMAKE_C_COMPILER`` and ``CMAKE_CXX_COMPILER`` need to be consistent with the AMReX
   make variable COMP to ensure matching OMP runtime libraries for use with the OpenMP NVector. 

#. ``CUDA_ARCH`` must be set to the appropriate value for the GPU being targeted

#. For more detailed instructions for installing SUNDIALS with different flags and versions see
   the `SUNDIALS documentation <https://computing.llnl.gov/projects/sundials/sundials-software>`_.

#. In the ``GNUmakefile`` for the application which uses the interface to SUNDIALS, add
   ``USE_SUNDIALS = TRUE`` and ``SUNDIALS_ROOT=${INSTALL_PREFIX}``. Note that one must define the
   ``SUNDIALS_LIB_DIR`` make variable to point to the location where the libraries are installed
   if they are not installed in the default location which is ``${INSTALL_PREFIX}/lib64``.

#. If the application uses the SUNDIALS CVODE time integrator package, then the variable
   ``USE_CVODE_LIBS = TRUE`` should also be added in the ``GNUmakefile`` for the application.
   If the application used the SUNDIALS ARKode time integrator package, then the variable
   ``USE_ARKODE_LIBS = TRUE`` should be added.

#. Fortran 2003 interfaces for the pgi compilers are currently not supported.


Note that SUNDIALS can also be installed via Spack:

   ::
      
      spack install sundials+cuda+f2003+openmp
  

SUNDIALS 5 Tutorials
--------------------------

AMReX provides in the ``amrex/Tutorials/SUNDIALS`` directory.
