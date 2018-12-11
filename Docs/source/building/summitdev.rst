Building WarpX for Summit-dev (OLCF)
====================================

For the `Summit-dev cluster
<https://www.olcf.ornl.gov/tag/summitdev/>`__ at OLCF,
use the following commands to download the source code, and switch to the
correct branch:

::

    mkdir warpx_directory
    cd warpx_directory

    git clone https://github.com/ECP-WarpX/WarpX.git
    cd WarpX
    git checkout dev
    cd ..

    git clone https://bitbucket.org/berkeleylab/picsar.git
    cd picsar
    git checkout gpu
    cd ..

    git clone https://github.com/AMReX-Codes/amrex.git
    cd amrex
    git checkout development
    cd ..


Then, use the following set of commands to compile:

::

    module load pgi
    module load cuda
    make -j 4 USE_GPU=TRUE COMP=pgi
