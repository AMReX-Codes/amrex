Running on specific platforms
=============================

Running on Cori KNL at NERSC
----------------------------

The batch script below can be used to run a WarpX simulation on 2 KNL nodes on
the supercomputer Cori at NERSC. Replace descriptions between chevrons ``<>``
by relevant values, for instance ``<job name>`` could be ``laserWakefield``.

.. literalinclude:: ../../../Tools/batchScripts/batch_cori.sh
   :language: bash

To run a simulation, copy the lines above to a file ``batch_cori.sh`` and
run
::

  sbatch batch_cori.sh

to submit the job.

For a 3D simulation with a few (1-4) particles per cell using FDTD Maxwell
solver on Cori KNL for a well load-balanced problem (in our case laser
wakefield acceleration simulation in a boosted frame in the quasi-linear
regime), the following set of parameters provided good performance:

* ``amr.max_grid_size=64`` and ``amr.blocking_factor=64`` so that the size of
  each grid is fixed to ``64**3`` (we are not using load-balancing here).

* **8 MPI ranks per KNL node**, with ``OMP_NUM_THREADS=8`` (that is 64 threads
  per KNL node, i.e. 1 thread per physical core, and 4 cores left to the
  system).

* **2 grids per MPI**, *i.e.*, 16 grids per KNL node.

Running on Summit at OLCF
-------------------------

The batch script below can be used to run a WarpX simulation on 2 nodes on
the supercomputer Summit at OLCF. Replace descriptions between chevrons ``<>``
by relevalt values, for instance ``<input file>`` could be
``plasma_mirror_inputs``. Note that the only option so far is to run with one
MPI rank per GPU.

.. literalinclude:: ../../../Tools/batchScripts/batch_summit.sh
   :language: bash

To run a simulation, copy the lines above to a file ``batch_summit.sh`` and
run
::

  bsub batch_summit.sh

to submit the job.

For a 3D simulation with a few (1-4) particles per cell using FDTD Maxwell
solver on Summit for a well load-balanced problem (in our case laser
wakefield acceleration simulation in a boosted frame in the quasi-linear
regime), the following set of parameters provided good performance:

* ``amr.max_grid_size=256`` and ``amr.blocking_factor=128``.

* **One MPI rank per GPU** (e.g., 6 MPI ranks for the 6 GPUs on each Summit
  node)

* **Two `128x128x128` grids per GPU**, or **one `128x128x256` grid per GPU**.

A batch script with more options regarding profiling on Summit can be found at
:download:`Summit batch script<../../../Tools/batchScripts/script_profiling_summit.sh>`
