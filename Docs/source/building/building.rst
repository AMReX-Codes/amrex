Building/installing WarpX
=========================

WarpX can be built with various options. This page describes the most basic
build, and points to instructions for more advanced builds.

Even if you are interested in more advanced builds, we recommend reading this
page first.

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

and to the branch ``dev`` of WarpX

::

    cd WarpX/
    git checkout dev
    cd ..

Basic compilation
-----------------

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

Advanced building instructions
------------------------------

.. toctree::
   :maxdepth: 1

   openpmd
   spectral
   gpu_local
   python
   spack

Building for specific plateforms
--------------------------------

.. toctree::
   :maxdepth: 1

   cori
   summit
