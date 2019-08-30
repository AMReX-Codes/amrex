Building WarpX for Summit (OLCF)
================================

For the `Summit cluster
<https://www.olcf.ornl.gov/summit/>`__ at OLCF,
use the following commands to download the source code, and switch to the
correct branch:

::

    mkdir warpx_directory
    cd warpx_directory

    git clone --branch dev https://github.com/ECP-WarpX/WarpX.git
    git clone --branch master https://bitbucket.org/berkeleylab/picsar.git
    git clone --branch development https://github.com/AMReX-Codes/amrex.git

Then, ``cd`` into the directory ``WarpX`` and use the following set of commands to compile:

::

    module load pgi
    module load cuda
    make -j 4 USE_GPU=TRUE COMP=pgi

See :doc:`../running_cpp/platforms` for more information on how to run
WarpX on Summit.
