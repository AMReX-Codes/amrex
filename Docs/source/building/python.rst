Installing WarpX as a Python package
------------------------------------

WarpX' Python bindings depend on ``numpy``, ``mpi4py``, and ``picmistandard``.

Type

.. code-block:: bash

   make -j 4 USE_PYTHON_MAIN=TRUE

or edit the ``GNUmakefile`` and set ``USE_PYTHON_MAIN=TRUE``, and type

.. code-block:: bash

   make -j 4

This will compile the code, and install the Python bindings as a package (named ``pywarpx``) in your standard Python installation (i.e. in your ``site-packages`` directory).
The note on compiler options from the previous section also holds when compiling the Python package.

In case you do not have write permissions to the default Python installation (e.g. typical on computer clusters), use the following command instead:

.. code-block:: bash

   make -j 4 PYINSTALLOPTIONS=--user

In this case, you can also set the variable ``PYTHONUSERBASE`` to set the folder where ``pywarpx`` will be installed.
