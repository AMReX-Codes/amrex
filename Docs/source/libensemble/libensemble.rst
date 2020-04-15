Run LibEnsemble on WarpX
========================

`LibEnsemble <https://github.com/Libensemble>`__ is a library to coordinate the concurrent evaluation of dynamic ensembles of calculations.
While a WarpX simulation can provide insight in some physics, it remains a single point evaluation in the space of parameters.
If you have a simulation ready for use, but would like to (i) scan over some input parameters uniformly for, e.g., a tolerance study, or (ii) have a random evaluation of the space of input parameters within a given span or (iii) tune some input parameters to optimize an output parameter, e.g., beam emittance, energy spread, etc., LibEnsemble provides these capabilities and will take care of tasks monitoring with fault tolerance on multiple platform (LibEnsemble targets modern HPC platforms like Summit).

Scripts to run LibEnsemble on WarpX simulations can be found in ``WarpX/Tools/LibEnsemble/``.
This documentation does not aim at giving a training on LibEnsemble, so please refer to the `LibEnsemble documentation <https://libensemble.readthedocs.io/en/develop/>`__ for technical details.

WarpX example problem for LibEnsemble study
-------------------------------------------

The WarpX example is built on a 2D input file (so that 1 simulation take < 1 min) of a 2-stage laser-wakefield simulation in a boosted frame.
It aims at optimizing emittance preservation in the coupling between two consecutive plasma accelerator stages.
Each stage accelerates an externally-injected electron beam, which charge has been tuned to show a decent amount of Ez field flattening due to longitudinal beam loading.
Each stage has a parabolic transverse profile to guide the laser pulse, and a uniform longitudinal profile with cos-shape ramps and the entrance and at the exit.
A fresh laser pulse is introduced at the entrance of each stage and deleted at the exit.
The beam transverse distribution is matched to the first stage and the beam is injected at the beginning of the plateau of the first stage, so that the emittance is conserved in the first stage.
The two stages are separated by a few-cm gap, and a focusing lens is placed in-between.
Note that this is a very low resolution simulation to show an example, so it is **not** close to numerical convergence.

In this example, we allow LibEnsemble to tune four **input parameters**:

  - Length of the downramp of the first stage

  - Longitudinal position of the focusing lens (between the two stages)

  - Strength of the focusing lens

  - Length of the downramp of the second stage

The **output parameter** that LibEnsemble minimizes is the beam emittance at the exit of the second stage, while making sure the charge loss is small.

The scripts provided can run on a local machine or on the Summit supercomputer at OLCF.
Two options are available: random sampling of parameter space or optimization on the output parameter.
For the latter, we are using the Asynchronously Parallel Optimization Solver for finding Multiple Minima `APOSMM <https://libensemble.readthedocs.io/en/develop/examples/gen_funcs.html#module-aposmm>`__ method provided by LibEnsemble.

Install LibEnsemble
-------------------

Besides a working WarpX executable, you have to install libEnsemble and its dependencies.

You can either install all packages via `conda` (recommended),

.. code-block:: sh

   conda install -c conda-forge libensemble matplotlib numpy scipy yt

or try to install the same dependencies via `pip` (pick one *or* the other):

.. literalinclude:: ../../../Tools/LibEnsemble/requirements.txt


What's in ``Tools/LibEnsemble``?
--------------------------------

See the `LibEnsemble User Guide <https://libensemble.readthedocs.io/en/develop/overview_usecases.html>`__ for an overview of LibEnsemble concepts.
In a nutshell, a user needs to define

  - A generator function ``gen_f`` that will generate inputs for the simulation, which can be done uniformly, randomly or using an optimizer.
    The generator output, i.e., the simulation input, is called ``'x'``.
    The generator is provided  by LibEnsemble.
    When the generator is an optimizer, it takes the simulation output called ``'f'`` as an input.

  - A simulation function ``sim_f`` that will take ``'x'`` as an input and return a single output value ``'f'``.
    In our example, ``sim_f`` modifies the WarpX input file depending on ``'x'``, launches a WarpX simulation and reads the simulation output plotfile to extract ``'f'``.

  - An allocator function ``alloc_f`` that will feed available workers with tasks.
    This is provided by LibEnsemble.

The files in ``Tools/LibEnsemble/`` are:

``run_libensemble_on_warpx.py``
   This is the main LibEnsemble script.
   It imports ``gen_f`` and ``alloc_f`` from LibEnsemble, ``sim_f`` from file ``warpx_simf.py`` (see below), defines dictionaries for parameters of each of these objects (``gen_specs`` includes lower and upper bound of each element in the input array ``'x'``, ``alloc_specs`` and ``sim_specs`` respectively) and runs LibEnsemble.

``warpx_simf.py``
   defines the ``sim_f`` function called ``run_warpx``:

   .. doxygenfunction:: run_warpx

``sim/inputs``
   WarpX input file. Some of its parameters are modified in ``run_warpx``.

``write_sim_input.py``
   (util) update one WarpX input file depending on values in ``'x'`` for this run.

   .. doxygenfunction:: write_sim_input

``read_sim_output.py``
   (util) Read WarpX plotfile and return ``'f'``.

   .. doxygenfunction:: read_sim_output

``plot_results.py``
   (util) Read LibEnsemble output files ``.npy`` and ``.pickle`` and plot output ``'f'`` (and other output, just for convenience) as a function of input from all simulations.

``all_machine_specs.py``
   (util) dictionaries of machine-specific parameters.
   For convenience, the maximum number of WarpX runs is defined there.

``summit_submit_mproc.sh``
   Submission script for LibEnsemble+WarpX on Summit.
   Make sure to edit this file and add your project ID for allocations.

Run the example
---------------

On Summit or for a local run, LibEnsemble can run with a Random Number Generator (easiest, good as a first try) or with an optimizer (requires python package ``nlopt``).
This is set by the variable ``generator_type`` in ``run_libensemble_on_warpx.py``.
We hereafter assume that all Python modules required are installed and that a WarpX executable is available.

Run locally
^^^^^^^^^^^

Adjust the ``local_specs`` dictionary in ``all_machine_specs.py`` to fix the path to the WarpX executable (and optionally change the number of cores and OpenMP threads), and run

.. code-block:: sh

    python run_libE_on_warpx.py --comms local --nworkers 3

This is adapted to a 4-core machine, as it will use:

- 1 process to run LibEnsemble

- 1 process (among the 3 workers) to run the generator

- 2 processes to run 2 concurrent simulations

Run on Summit at OLCF
^^^^^^^^^^^^^^^^^^^^^

.. code-block:: sh

    bsub summit_submit_mproc.sh
