.. _building-openpmd:

Building WarpX with support for openPMD output
==============================================

WarpX can dump data in the `openPMD format <https://github.com/openPMD>`__.
This feature currently requires to have a parallel version of HDF5 installed ;
therefore we recommend to use `spack <https://
spack.io>`__ in order to facilitate the installation.

More specifically, we recommend that you try installing the
`openPMD-api library 0.11.0a or newer <https://openpmd-api.readthedocs.io/en/0.11.0-alpha/>`__
using spack (first section below). If this fails, a back-up solution
is to install parallel HDF5 with spack, and then install the openPMD-api
library from source.

In order to install spack, you can simply do:

.. code-block:: bash

   git clone https://github.com/spack/spack.git
   export SPACK_ROOT=$PWD/spack
   . $SPACK_ROOT/share/spack/setup-env.sh

You may want to auto-activate spack when you open a new terminal by adding this to your ``$HOME/.bashrc`` file:

.. code-block:: bash

   echo -e "# activate spack package manager\n. ${SPACK_ROOT}/share/spack/setup-env.sh" >> $HOME/.bashrc


WarpX Development Environment with Spack
----------------------------------------

Create and activate a Spack environment with all software needed to build WarpX

.. code-block:: bash

   spack env create warpx-dev    # you do this once
   spack env activate warpx-dev
   spack add gmake
   spack add mpi
   spack add openpmd-api
   spack add pkg-config
   spack install

This will download and compile all dependencies.

Whenever you need this development environment in the future, just repeat the quick ``spack env activate warpx-dev`` step.
For example, we can now compile WarpX by ``cd``-ing into the ``WarpX`` folder and typing:

.. code-block:: bash

   spack env activate warpx-dev
   make -j 4 USE_OPENPMD=TRUE

You will also need to load the same spack environment when running WarpX, for instance:

.. code-block:: bash

   spack env activate warpx-dev
   mpirun -np 4 ./warpx.exe inputs

You can check which Spack environments exist and if one is still active with

.. code-block:: bash

   spack env list  # already created environments
   spack env st    # is an environment active?


Installing openPMD-api from source
----------------------------------

You can also build openPMD-api from source, e.g. to build against the module environment of a supercomputer cluster.

First, load the according modules of the cluster to support the openPMD-api dependencies.
You can find the `required and optional dependencies here <https://github.com/openPMD/openPMD-api#dependencies_`.

You usually just need a C++ compiler, CMake, and one or more file backend libraries, such as HDF5 and/or ADIOS2.
See for example `our installation guidelines for Cori :ref`<building-cori-openPMD>`.

Then, in the ``$HOME/warpx_directory/``, download and build openPMD-api:

.. code-block:: bash

   git clone https://github.com/openPMD/openPMD-api.git
   mkdir openPMD-api-build
   cd openPMD-api-build
   cmake ../openPMD-api -DopenPMD_USE_PYTHON=OFF -DCMAKE_INSTALL_PREFIX=$HOME/warpx_directory/openPMD-install/ -DCMAKE_INSTALL_RPATH_USE_LINK_PATH=ON -DCMAKE_INSTALL_RPATH='$ORIGIN'
   cmake --build . --target install

Finally, compile WarpX:

.. code-block:: bash

   cd ../WarpX
   export PKG_CONFIG_PATH=$HOME/warpx_directory/openPMD-install/lib/pkgconfig:$PKG_CONFIG_PATH
   export CMAKE_PREFIX_PATH=$HOME/warpx_directory/openPMD-install:$CMAKE_PREFIX_PATH

   make -j 4 USE_OPENPMD=TRUE

When running WarpX, we will recall where you installed openPMD-api via RPATHs, so you just need to load the same module environment as used for building (same MPI, HDF5, ADIOS2, for instance).

.. code-block:: bash

   # module load ...  (compiler, MPI, HDF5, ADIOS2, ...)

   mpirun -np 4 ./warpx.exe inputs
