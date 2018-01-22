.. role:: cpp(code)
   :language: c++


Downloading the Code
====================

The source code of  is available at
https://github.com/AMReX-Codes/amrex. The GitHub repo is our
central repo for development. The development branch
includes the latest state of the code, and it is merged into the
master branch on a monthly basis. The master branch is
considered the release branch. The releases are tagged with version
number YY.MM (e.g., 17.04). The MM part of the
version is incremented every month, and the YY part every year.
Bug fix releases are tagged with YY.MM.patch (e.g.,
17.04.1).

Example: Hello World
====================

The source code of this example is at
amrex/Tutorials/Basic/HelloWorld_C/ and is also shown below.

.. highlight:: c++

::

     #include <AMReX.H>
     #include <AMReX_Print.H>

     int main(int argc, char* argv[])
     {
         amrex::Initialize(argc,argv);
         amrex::Print() << "Hello world from AMReX version " 
                        << amrex::Version() << "\n";
         amrex::Finalize();
     }

The main body of this short example contains three statements.
Usually the first and last statements for the :cpp:`int main(...)` function of
every program should be calling :cpp:`amrex::Initialize` and :cpp:`Finalize`, 
respectively. The second statement calls `amrex::Print`` to print out
a string that includes the AMReX version returned by the :cpp:`amrex::Version`
function. The example code includes two AMReX header files. Note
that the name of all AMReX header files starts with ``AMReX_``
(or just AMReX in the case of AMReX.H). All AMReX C++ functions are in the 
:cpp:`amrex` namespace.

Building the Code
-----------------

You build the code in the amrex/Tutorials/Basic/HelloWorld_C/
directory. Typing ``make`` will start the compilation process and
result in an executable named ``main3d.gnu.DEBUG.ex``. The name
shows that the GNU compiler with debug options set by AMReX is used.
It also shows that the executable is built for 3D. Although this
simple example code is dimension independent, dimensionality does matter
for all non-trivial examples. The build process can be adjusted by
modifying the amrex/Tutorials/Basic/HelloWorld_C/GNUmakefile file.
More details on how to build AMReX can be found in :ref:`Chap:BuildingAMReX`.

Running the Code
----------------

The example code can be run as follows,

.. highlight:: console

::

      ./main3d.gnu.DEBUG.ex

The result may look like,

.. highlight:: console

::

      Hello world from AMReX version 17.05-30-g5775aed933c4-dirty

The version string means the current commit 5775aed933c4 (note
that the first letter g in g577.. is not part of the hash)
is based on 17.05 with 30 additional commits and the AMReX work tree 
is dirty (i.e. there are uncommitted changes).

In the GNUmakefile there are compilation options for DEBUG
mode (less optimized code with more error checking), dimensionality,
compiler type, and flags to enable MPI and/or OpenMP parallelism.
If there are multiple instances of a parameter, the last instance
takes precedence.

Parallelization
---------------

Now let's build with MPI by typing ``make USE_MPI=TRUE`` (alternatively
you can set ``USE_MPI=TRUE`` in the GNUmakefile). This
should make an executable named ``main3d.gnu.DEBUG.MPI.ex``. Note
MPI in the file name. You can then run,

.. highlight:: console

::

      mpiexec -n 4 ./main3d.gnu.DEBUG.MPI.ex

The result may look like,

.. highlight:: console

::

      MPI initialized with 4 MPI processes
      Hello world from AMReX version 17.05-30-g5775aed933c4-dirty

If the compilation fails, you are referred to :ref:`Chap:BuildingAMReX` 
for more details on how to configure the build system.

If you want to build with OpenMP, type make ``USE_OMP=TRUE``.
This should make an executable named ``main3d.gnu.DEBUG.OMP.ex``. Note
OMP in the file name. Make sure the ``OMP_NUM_THREADS``
environment variable is set on your system. You can then run,

.. highlight:: console

::

      ./main3d.gnu.DEBUG.OMP.ex

The result may look like,

.. highlight:: console

::

      OMP initialized with 4 OMP threads
      Hello world from AMReX version 17.06-287-g51875485fe51-dirty

Note that you can build with both USE_MPI=TRUE and USE_OMP=TRUE.
You can then run,

.. highlight:: console

::

      mpiexec -n 2 ./main3d.gnu.DEBUG.MPI.OMP.ex

The result may look like,

.. highlight:: console

::

      MPI initialized with 2 MPI processes
      OMP initialized with 4 OMP threads
      Hello world from AMReX version 17.06-287-g51875485fe51-dirty

.. _sec:heat equation:

Example: Heat Equation Solver
=============================

We now look at a more complicated example at
amrex/Tutorials/Basic/HeatEquation_EX1_C and show how simulation
results can be visualized. This example solves the heat equation,

.. math:: \frac{\partial\phi}{\partial t} = \nabla^2\phi

using forward Euler temporal integration on a periodic domain.
We could use a 5-point (in 2D) or 7-point (in 3D) stencil, but for demonstration
purposes we spatially discretize the PDE by first constructing fluxes on cell faces, e.g.,

.. math:: F_{i+^1\!/_2,\,j} = \frac{\phi_{i+1,j}-\phi_{i,j}}{\Delta x},

and then taking the divergence to update the cells,

.. math::

   \phi_{i,\,j}^{n+1} = \phi_{i,\,j}^n 
   + \frac{\Delta t}{\Delta x}\left(F_{i+^1\!/_2,\,j}-F_{i-^1\!/_2,\,j}\right)
   + \frac{\Delta t}{\Delta y}\left(F_{i,\,j+^1\!/_2}-F_{i,\,j-^1\!/_2}\right)

Don't worry about the implementation details of the code.
You will be able to understand the code in this example after
reading the section on the :ref:`Chap:Basics`.

Building and Running the Code
-----------------------------

To build a 2D executable, type ``make DIM=2``. This will generate
an executable named ``main2d.gnu.ex``. To run it, type,

.. highlight:: console

::

      ./main2d.gnu.DEBUG.ex inputs_2d

Note that the command takes a file inputs_2d. When the run
finishes, you will have a number of plotfiles, plt00000,
plt01000, etc. The calculation solves the heat equation in 2D on a
:math:`256 \times 256` cells domain. It runs :math:`10,000` steps and makes a
plotfile every :math:`1,000` steps. These are runtime parameters that can
be adjusted in inputs_2d.

Visualization
=============

There are several visualization tools that can be used for  plotfiles. The standard tool used within the
AMReX-community is AmrVis, a package developed and supported
by CCSE that is designed specifically for highly efficient visualization
of block-structured hierarchical AMR data.
Plotfiles can also be viewed using the VisIt, ParaView, and yt packages.
Particle data can be viewed using ParaView.
Refer to Chapter on :ref:`Chap:Visualization` for how to use each of these tools.
