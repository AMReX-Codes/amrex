Building WarpX for Summit (OLCF)
================================

For the `Summit cluster
<https://www.olcf.ornl.gov/summit/>`__ at OLCF,
use the following commands to download the source code, and switch to the
correct branch:

::

    mkdir warpx_directory
    cd warpx_directory

    git clone https://github.com/ECP-WarpX/WarpX.git
    git clone https://bitbucket.org/berkeleylab/picsar.git
    git clone --branch development https://github.com/AMReX-Codes/amrex.git

Then, ``cd`` into the directory ``WarpX`` and use the following set of commands to compile:

::

    module load gcc
    module load cuda
    make -j 4 USE_GPU=TRUE

See :doc:`../running_cpp/platforms` for more information on how to run WarpX on Summit.

See :doc:`../visualization/yt` for more information on how to visualize the simulation results. In order to build yt **and** a Python environment, you can download the yt installation script as explained in section "All-in-One Installation Script" of the `yt installation web page <https://yt-project.org/doc/installing.html>`__, and run

.. code-block:: sh

    module purge
    module load gcc

Then modify the few first lines of the installation script ``install_script.sh`` to have

::

   INST_YT_SOURCE=1
   INST_SCIPY=1

and follow the instructions from the yt website.
