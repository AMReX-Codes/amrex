Building/installing WarpX
=========================

WarpX can currently be built (and run) in two variants:
   - as a compiled executable (run with the command line)
   - as a Python package (run through a Python script)

Currently, for both of these options, the user needs to build the code from source.

Downloading the source code
---------------------------

Clone the source codes of WarpX, and its dependencies AMReX and PICSAR into one
single directory (e.g. ``warpx_directory``):

::

    mkdir warpx_directory
    cd warpx_directory
    git clone https://github.com/ECP-WarpX/WarpX.git
    git clone https://bitbucket.org/berkeleylab/picsar.git
    git clone https://github.com/AMReX-Codes/amrex.git

Then switch to the branch ``development`` of AMReX

::

    cd amrex/
    git checkout development
    cd ..

Compiling WarpX as an executable
--------------------------------

``cd`` into the directory ``WarpX`` and type

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


Installing WarpX as a Python package
------------------------------------

Type

::

    make -j 4 USE_PYTHON_MAIN=TRUE

or edit the GNUmakefile and set `USE_PYTHON_MAIN=TRUE`, and type

::

    make -j 4

This will compile the code, and install the Python bindings as a package (named
``pywarpx``) in your standard Python installation (i.e. in your
``site-packages`` directory). The note on compiler options from the previous
section also holds when compiling the Python package.

In case you do not have write permissions to the default Python installation (e.g. typical on computer clusters), use the following command instead:

::

   make -j 4 PYINSTALLOPTIONS=--user

In this case, you can also set the variable `PYTHONUSERBASE` to set the folder where `pywarpx` will be installed.

Advanced building instructions
------------------------------


.. toctree::
   :maxdepth: 1

   spectral
   cori
   summitdev
