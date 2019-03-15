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

    * with the GNU compiler

    ::

        module swap PrgEnv-intel PrgEnv-gnu
        make -j 16 COMP=gnu

In order to compile for the **Knight's Landing (KNL) architecture**:

    * with the Intel compiler

    ::

        module swap craype-haswell craype-mic-knl
        make -j 16 COMP=intel

    * with the GNU compiler

    ::

        module swap craype-haswell craype-mic-knl
        module swap PrgEnv-intel PrgEnv-gnu
        make -j 16 COMP=gnu


Building WarpX with openPMD support
-----------------------------------

First, load the appropriate modules:

::

    module swap craype-haswell craype-mic-knl
    module swap PrgEnv-intel PrgEnv-gnu
    module load cmake/3.11.4
    module load cray-hdf5-parallel

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
    export OPENPMD_LIB_PATH=../openPMD-install/lib64
    export OPENPMD_INCLUDE_PATH=../openPMD-install/include
    make -j 16 COMP=gnu USE_OPENPMD=TRUE
