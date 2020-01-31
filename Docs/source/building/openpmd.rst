.. _building-openpmd:

Building WarpX with support for openPMD output
==============================================

WarpX can dump data in the `openPMD format <https://github.com/openPMD>`__.
This feature currently requires to have a parallel version of HDF5 installed ;
therefore we recommend to use `spack <https://
spack.io>`__ in order to facilitate the installation.

More specifically, we recommend that you try installing the
`openPMD-api library 0.10.3a or newer <https://openpmd-api.readthedocs.io/en/0.10.3-alpha/>`__
using spack (first section below). If this fails, a back-up solution
is to install parallel HDF5 with spack, and then install the openPMD-api
library from source.

In order to install spack, you can simply do:

::

  git clone https://github.com/spack/spack.git
  export SPACK_ROOT=/path/to/spack
  . $SPACK_ROOT/share/spack/setup-env.sh

(You may want to add the last 2 lines to your ``.bashrc`` file.)


Building openPMD support, by installing openPMD-api directly from spack
-----------------------------------------------------------------------

First, install the openPMD-api library:

::

    spack install openpmd-api -python

Then, ``cd`` into the ``WarpX`` folder, and type:

::

    spack load -r openpmd-api
    make -j 4 USE_OPENPMD=TRUE

You will also need to load the same spack environment when running WarpX, for instance:

::

    spack load -r openpmd-api

    mpirun -np 4 ./warpx.exe inputs

Building openPMD support, by installing openPMD-api from source
---------------------------------------------------------------

First, install the openPMD-api library, and load it in your environment:

::

    spack install hdf5
    spack install adios2
    spack load -r hdf5
    spack load -r adios2

Then, in the `warpx_directory`, download and build the openPMD API:

::

    git clone https://github.com/openPMD/openPMD-api.git
    mkdir openPMD-api-build
    cd openPMD-api-build
    cmake ../openPMD-api -DopenPMD_USE_PYTHON=OFF -DCMAKE_INSTALL_PREFIX=../openPMD-install/ -DCMAKE_INSTALL_RPATH_USE_LINK_PATH=ON -DCMAKE_INSTALL_RPATH='$ORIGIN'
    cmake --build . --target install

Finally, compile WarpX:

::

    cd ../WarpX
    export PKG_CONFIG_PATH=$PWD/../openPMD-install/lib/pkgconfig:$PKG_CONFIG_PATH
    make -j 4 USE_OPENPMD=TRUE

You will also need to load the same spack environment when running WarpX, for instance:

::

    spack load -r hdf5
    spack load -r adios2

    mpirun -np 4 ./warpx.exe inputs
