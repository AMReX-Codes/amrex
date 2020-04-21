Using WarpX with Python
===========================

WarpX uses the PICMI standard for its Python input files.
See `PICMI information and source code <https://github.com/picmi-standard/picmi>`__.

WarpX can be run in one of two modes. It can run as a preprocessor, using the
Python input file to generate an input file to be used by the C++ version, or
it can be run directly from Python.

Example input files can be found in the section :doc:`examples`.
The examples support running in both modes by commenting and uncommenting the appropriate lines.


Using Python input as a preprocessor
------------------------------------

In this case, only the pure Python version needs to be installed, as described in :doc:`../building/python.rst`.

In order to run a new simulation:

* Create a **new directory**, where the simulation will be run.

* Add a **Python script** in the directory.

The input file should have the line like "sim.write_input_file(file_name = 'inputs_from_PICMI')"
which runs the preprocessor, generating the C++ input file.

* **Run** the script with Python:

::

    python <python_script>

where ``<python_script>`` is the name of the script.
This creates the WarpX input file that you can run as normal with the WarpX executable.


Running WarpX directly from Python
----------------------------------

For this, a full Python installation of WarpX is required, as described in :doc:`../building/python.rst`.

In order to run a new simulation:

* Create a **new directory**, where the simulation will be run.

* Add a **Python script** in the directory.

The input file should have the line "sim.step()" which runs the simulation.

* **Run** the script with Python:

::

    mpirun -np <n_ranks> python <python_script>

where ``<n_ranks>`` is the number of MPI ranks used, and ``<python_script>``
is the name of the script.

