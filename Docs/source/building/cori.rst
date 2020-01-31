Building WarpX for Cori (NERSC)
===============================

Standard build
--------------

For the `Cori cluster
<http://www.nersc.gov/users/computational-systems/cori/>`__ at NERSC,
you need to type the following command when compiling:

.. note::

   In order to compile the code with a spectral solver, type

   ::

       module load cray-fftw

   before typing any of the commands below, and add ``USE_PSATD=TRUE``
   at the end of the command containing ``make``.

In order to compile for the **Haswell architecture**:

    * with the Intel compiler

    ::

        make -j 16 COMP=intel

    * with the GNU compiler (``COMP=gcc`` is the corresponding optional flag)

    ::

        module swap PrgEnv-intel PrgEnv-gnu
        make -j 16 COMP=gcc

In order to compile for the **Knight's Landing (KNL) architecture**:

    * with the Intel compiler

    ::

        module swap craype-haswell craype-mic-knl
        make -j 16 COMP=intel

    * with the GNU compiler

    ::

        module swap craype-haswell craype-mic-knl
        module swap PrgEnv-intel PrgEnv-gnu
        make -j 16 COMP=gcc

See :doc:`../running_cpp/platforms` for more information on how to run
WarpX on Cori.

GPU Build
---------

To compile on the experimental GPU nodes on Cori, you first need to purge
your modules, most of which won't work on the GPU nodes.

::

    module purge

Then, you need to load the following modules:

::

    module load modules esslurm gcc/7.3.0 cuda mvapich2

You can also use OpenMPI-UCX instead of mvapich: openmpi/4.0.1-ucx-1.6.

Then, you need to use slurm to request access to a GPU node:

::

    salloc -C gpu -N 1 -t 30 -c 10 --gres=gpu:1 -A m1759

This reserves 10 logical cores (5 physical), 1 GPU.
The latest documentation can be found here: https://docs-dev.nersc.gov/cgpu/access
Note that you can't cross-compile for the GPU nodes - you have to log on to one
and then build your software.

Finally, navigate to the base of the WarpX repository and compile in GPU mode:

::

    make -j 16 USE_GPU=TRUE


Building WarpX with openPMD support
-----------------------------------

First, load the appropriate modules:

::

    module swap craype-haswell craype-mic-knl
    module swap PrgEnv-intel PrgEnv-gnu
    module load cmake/3.14.4
    module load cray-hdf5-parallel
    module load adios/1.13.1
    export CRAYPE_LINK_TYPE=dynamic

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
    export PKG_CONFIG_PATH=$PWD/../openPMD-install/lib64/pkgconfig:$PKG_CONFIG_PATH
    make -j 16 COMP=gcc USE_OPENPMD=TRUE

In order to run WarpX, load the same modules again.
