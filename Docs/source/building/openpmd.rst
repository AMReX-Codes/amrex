Building WarpX with support for openPMD output
==============================================

WarpX can dump data in the `openPMD format <https://github.com/openPMD>`__.
This feature currently requires to have a parallel version of HDF5 installed ;
therefore we recommend to use `spack <https://
spack.io>`__ in order to facilitate the installation.

More specifically, we recommend that you try installing the
`openPMD-api library <https://openpmd-api.readthedocs.io/en/0.8.0-alpha/>`__
using spack (first section below). If this fails, a back-up solution
is to install parallel HDF5 with spack, and then install the openPMD-api
library from source.

Building openPMD support, by installing openPMD-api directly from spack
-----------------------------------------------------------------------

First, install the openPMD-api library:

::

    spack install openpmd-api -shared -json -python ^hdf5+mpi ^openmpi

Then, ``cd`` into the ``WarpX`` folder, and type:

::

    spack load openmpi
    spack load hdf5
    spack load openpmd-api
    make -j 4 USE_OPENPMD=TRUE

You will also need to load the same spack environment when running WarpX, for instance:

::

    spack load openmpi
    spack load hdf5
    spack load openpmd-api

    mpirun -np 4 ./warpx.exe inputs

Building openPMD support, by installing openPMD-api from source
---------------------------------------------------------------

First, install the openPMD-api library, and load it in your environment:

::

    spack install hdf5+mpi ^openmpi
    spack load openmpi
    spack load hdf5

Then, in the `warpx_directory`, download and build the openPMD API:

::

    git clone https://github.com/openPMD/openPMD-api.git
    mkdir openPMD-api-build
    cd openPMD-api-build
    cmake ../openPMD-api -DopenPMD_USE_PYTHON=OFF -DopenPMD_USE_JSON=OFF -DCMAKE_INSTALL_PREFIX=../openPMD-install/ -DBUILD_SHARED_LIBS=OFF
    cmake --build . --target install

Finally, compile WarpX:

::

    cd ../WarpX
    export OPENPMD_LIB_PATH=../openPMD-install/lib
    export OPENPMD_INCLUDE_PATH=../openPMD-install/include
    make -j 4 USE_OPENPMD=TRUE

You will also need to load the same spack environment when running WarpX, for instance:

::

    spack load openmpi
    spack load hdf5

    mpirun -np 4 ./warpx.exe inputs
