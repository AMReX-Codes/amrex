Building WarpX with support for openPMD output
==============================================

WarpX can dump data in the `openPMD format <https://github.com/openPMD>`__.
This feature currently requires to have a parallel version of HDF5 installed ;
therefore we recommend to use `spack <https://
spack.io>`__ in order to facilitate the installation.

In order to build WarpX with openPMD support, you first need to install
the openPMD-api library:

::

    spack spec openpmd-api -shared -json -python ^hdf5+mpi ^openmpi

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
