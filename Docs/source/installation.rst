Installation
============

Downloading the code
~~~~~~~~~~~~~~~~~~~~

You should clone the source codes of AMReX, PICSAR and WarpX into one
single directory (e.g. ``warpx_directory``):

::

    mkdir warpx_directory
    cd warpx_directory
    git clone https://bitbucket.org/berkeleylab/amrex.git
    git clone https://bitbucket.org/berkeleylab/picsar.git
    git clone https://bitbucket.org/berkeleylab/warpx.git

You should then switch to the branch ``development`` of AMReX

::

    cd amrex/
    git checkout development
    cd ..

Compiling the code
~~~~~~~~~~~~~~~~~~

``cd`` into the directory ``warpx`` and type

::

    make -j 4

(in order to compile the code in parallel on 4 cores).  This will
generate an executable file in the ``Bin`` directory.

In order to clean a previously compiled version:

::

    make realclean

Running the Langmuir tests
--------------------------

The folder ``tests/Langmuir`` contains code that allow the user to run a
Langmuir wave case with either WarpX or PICSAR, and to compare the
results. The instructions below explain how to do this.
