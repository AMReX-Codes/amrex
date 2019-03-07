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
    cd WarpX
    git checkout dev
    cd ..

    git clone https://bitbucket.org/berkeleylab/picsar.git
    cd picsar
    git checkout master
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

In order to submit a simulation, create a file `submission_script` with
the following text (replace bracketed variables):

::

    #!/bin/bash
    #BSUB -J <jobName>
    #BSUB -W <requestedTime>
    #BSUB -nnodes <numberOfNodes>
    #BSUB -P <accountNumber>

    module load pgi
    module load cuda

    omp=1
    export OMP_NUM_THREADS=${omp}
    jsrun -n <numberOfNodes> -a 6 -g 6 -c 6 --bind=packed:${omp} --smpiargs="-gpu" <warpxExecutable> <inputScript>


Then run

::

    bsub submission_script
