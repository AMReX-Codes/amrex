Running WarpX as a compiled executable
======================================

After compiling the code (see the section :doc:`../installation`), the
WarpX executable is stored in the folder ``warpx/Bin``. (Its name starts
with `main` but depends on the compiler options.)

How to run a new simulation
~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Create a **new directory**, where the simulation will be run.
* Copy the **executable** to this directory:

::

    cp warpx/Bin/<warpx_executable> <run_directory>/warpx.exe


where ``<warpx_executable>`` should be replaced by the actual name
of the executable (see above) and ``<run_directory>`` by the actual
path to the run directory.


* Add an **input file** in the directory.

This file contains the numerical and physical parameters that define
the situation to be simulated.
Example input files can be found in the section :doc:`examples`.
The different parameters in these files are explained in the section
:doc:`parameters`.


* **Run** the executable:

::

    mpirun -np <n_ranks> ./warpx.exe <input_file>

where ``<n_ranks>`` is the number of MPI ranks used, and ``<input_file>``
is the name of the input file.


For more details
~~~~~~~~~~~~~~~~

.. toctree::
   :maxdepth: 1

   examples
   parameters
