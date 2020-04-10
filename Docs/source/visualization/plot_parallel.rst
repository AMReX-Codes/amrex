Out-of-the-box plotting script
==============================

A ready-to-use python script for plotting simulation results is available at
:download:`plot_parallel.py<../../../Tools/PostProcessing/plot_parallel.py>`. Feel free to
use it out-of-the-box or to modify it to suit your needs.

Dependencies
------------

Most of its dependencies are standard Python packages, that come with a default
Anaconda installation or can be installed with ``pip`` or ``conda``:
`os, matplotlib, sys, argparse, matplotlib, scipy`.

Additional dependencies are ``yt >= 3.5`` ( or ``yt >= 3.6`` if you are using
rigid injection, see section :doc:`yt` on how to install ``yt``), and ``mpi4py``.

Run serial
----------

Executing the script with

::

    python plot_parallel.py

will loop through plotfiles named ``plt?????`` (e.g., ``plt00000``, ``plt00100`` etc.)
and save one image per plotfile. For a 2D simulation, a 2D colormap of the Ez
field is plotted by default, with 1/20 of particles of each species (with different colors).
For a 3D simulation, a 2D colormap of the central slices in `y` is plotted, and particles
are handled the same way.

The script reads command-line options (which field and particle species, rendering with
`yt` or `matplotlib`, etc.). For the full list of options, run

::

    python plot_parallel.py --help

In particular, option ``--plot_Ey_max_evolution`` shows you how to plot the evolution of
a scalar quantity over time (by default, the max of the Ey field). Feel free to modify it
to plot the evolution of other quantities.

Run parallel
------------

To execute the script in parallel, you can run for instance

::

    mpirun -np 4 python plot_parallel.py --parallel

In this case, MPI ranks will share the plotfiles to process as evenly as possible.
Note that each plotfile is still processed in serial. When option
``--plot_Ey_max_evolution`` is on, the scalar quantity is gathered to rank 0, and
rank 0 plots the image.

If all dependencies are satisfied, the script can be used on Summit or Cori. For
instance, the following batch script illustrates how to submit a post-processing
batch job on Cori haswell with some options:

.. literalinclude:: ../../../Tools/PostProcessing/cori_postproc_script.sh
    :language: bash
