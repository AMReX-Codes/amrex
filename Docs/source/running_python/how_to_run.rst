How to run a new simulation
===========================

After installing WarpX as a Python package (see the section
:doc:`../installation`), you can use its functionalities in a Python script
to run a simulation.

In order to run a new simulation:

* Create a **new directory**, where the simulation will be run.

* Add a **Python script** in the directory.

This file contains the numerical and physical parameters that define
the situation to be simulated.
Example input files can be found in the section :doc:`examples`.

* **Run** the script with Python:

::

    mpirun -np <n_ranks> python <python_script>

where ``<n_ranks>`` is the number of MPI ranks used, and ``<python_script>``
is the name of the script.
