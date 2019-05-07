Building WarpX with GPU support (Linux only)
--------------------------------------------

.. warning::

  In order to build WarpX on a specific GPU cluster (e.g. Summit),
  look for the corresponding specific instructions, instead
  of those on this page.

In order to build WarpX with GPU support, it necessary to install:

- The PGI compiler. The free PGI Community Edition can be installed from
  `this page <https://www.pgroup.com/products/community.htm>`__.

- `Spack <https://spack.readthedocs.io/en/latest/index.html>`__:

  ::

    git clone https://github.com/spack/spack.git
    export SPACK_ROOT=/path/to/spack
    . $SPACK_ROOT/share/spack/setup-env.sh

(You may want to add the last 2 lines to your ``.bashrc`` file.)

Next, you will to install of version of MPI which is compatible with the PGI
compiler (using Spack). First of all, make sure that Spack is aware of the PGI
compiler:

  ::

    spack compiler find
    spack compilers

Then run:

  ::

    spack install mvapich2 fabrics=sock %pgi ^libpciaccess%gcc

Finally, WarpX can be compiled with

  ::

    spack load mvapich2%pgi
    make -j 4 USE_GPU=TRUE COMP=pgi
