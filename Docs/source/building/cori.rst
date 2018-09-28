Building WarpX for Cori (NERSC)
===============================

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
