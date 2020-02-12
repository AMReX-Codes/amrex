Building WarpX with Spack
===============================

WarpX can be installed using Spack. From the Spack web page: "Spack is a package management tool designed to support multiple
versions and configurations of software on a wide variety of platforms and environments."

Spack is available from `github <https://github.com/spack/spack>`__. Spack only needs to be cloned and can be used right away - there are no installation
steps. The spack command, "spack/bin/spack", can be used directly or "spack/bin" can be added to your execute path.

WarpX is built with the single command

::

    spack install warpx

This will build the 3-D version of WarpX using the master branch.
At the very end of the output from build sequence, Spack tells you where the WarpX executable has been placed.
Alternatively, the "spack load" command can be configured so that "spack load warpx" will put the executable in your execute path.

Other variants of WarpX can be installed, for example

::

    spack install warpx dims=2

will build the 2-D version.

::

    spack install warpx debug=True

will build with debugging turned on.

::

    spack install warpx %intel

will build using the intel compiler (instead of gcc).

The Python verson of WarpX is not yet available with Spack.
