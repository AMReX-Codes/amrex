Installation
============

Downloading the code
~~~~~~~~~~~~~~~~~~~~

You should clone the source codes of AMReX, PICSAR and WarpX into one
single directory (e.g. ``warpx_directory``):

::

    mkdir warpx_directory
    cd warpx_directory
    git clone https://github.com/AMReX-Codes/amrex.git
    git clone https://bitbucket.org/berkeleylab/picsar.git
    git clone https://github.com/AMReX-Codes/amrex.git

You should then switch to the branch ``development`` of AMReX

::

    cd amrex/
    git checkout development
    cd ..

Compiling WarpX as an executable
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``cd`` into the directory ``warpx`` and type

::

    make -j 4

This will generate an executable file in the ``Bin`` directory.

.. note::

    The compilation options are set in the file ``GNUmakefile``. The default
    options correspond to an optimized code for 3D geometry. You can modify the
    options in this file in order to (for instance):

        * Use 2D geometry
        * Disable OpenMP
        * Profile or debug the code
        * Choose a given compiler

    For a description of these different options, see the `corresponding page
    <https://amrex-codes.github.io/amrex/BuildingAMReX.html#building-with-gnu-make>`__ in
    the AMReX documentation.

    Alternatively, instead of modifying the file ``GNUmakefile``, you can
    directly pass the options in command line ; for instance:

    ::

        make -j 4 USE_OMP=FALSE


In order to clean a previously compiled version:

::

    make realclean

In order to learn how to use the executable, see the section :doc:`running_cpp/running_cpp`.


Compiling WarpX as a Python package
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``cd`` into the directory ``warpx/Python`` (instead of ``warpx``
in the previous section) and type

::

    make -j 4

This will compile the code, and install the Python bindings as a package (named
``pywarpx``) in your standard Python installation (i.e. in your
``site-packages`` directory). The note on compiler options from the previous
section also hold when compiling the Python package.

In order to learn how to use the Python package, see the section :doc:`running_python/running_python`.


Compiling WarpX for Cori
~~~~~~~~~~~~~~~~~~~~~~~~

For the `Cori cluster
<http://www.nersc.gov/users/computational-systems/cori/>`__ at NERSC,
you need to type the following command when compiling:

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
