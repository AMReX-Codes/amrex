Building WarpX with GPU support (Linux only)
--------------------------------------------

.. warning::

  In order to build WarpX on a specific GPU cluster (e.g. Summit),
  look for the corresponding specific instructions, instead
  of those on this page.

In order to build WarpX with GPU support, make sure that you have `cuda`
and `mpich` installed on your system. (Compiling with `openmpi` currently
fails.) Then compile WarpX with the option `USE_GPU=TRUE`, e.g.

::

  make -j 4 USE_GPU=TRUE
