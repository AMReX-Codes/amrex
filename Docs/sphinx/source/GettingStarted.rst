In this chapter, we will walk you through two simple examples. It is
assumed here that your machine has GNU Make, Python, GCC (including
gfortran), and MPI, although  can be built with CMake and other
compilers.

Downloading the Code
====================

The source code of  is available at
https://github.com/AMReX-Codes/amrex. The GitHub repo is our
central repo for development. The development branch
includes the latest state of the code, and it is merged into the master branch on a monthly basis. The master branch is
considered the release branch. The releases are tagged with version
number YY.MM (e.g., 17.04). The MM part of the
version is incremented every month, and the YY part every year.
Bug fix releases are tagged with YY.MM.patch (e.g., 17.04.1).

Example: Hello World
====================

The source code of this example is at amrex/Tutorials/Basic/HelloWorld\_C/ and is also shown below.

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
Usually the first and last statements for the main function of
every program should be calling amrex:: and
, respectively. The second statement calls amrex:: to print out a string that includes the
 version returned by the amrex::
function. The example code includes two  header files. Note
that the name of all  header files starts with AMReX\_
(or just AMReX in the case of AMReX.H). All  functions are in the amrex namespace.

Building the Code
-----------------

You build the code in the amrex/Tutorials/Basic/HelloWorld\_C/
directory. Typing make will start the compilation process and
result in an executable named main3d.gnu.DEBUG.ex. The name
shows that the GNU compiler with debug options set by  is used.
It also shows that the executable is built for 3D. Although this
simple example code is dimension independent, the dimension matters
for all non-trivial examples. The build process can be adjusted by
modifying the amrex/Tutorials/Basic/HelloWorld\_C/GNUmakefile file.
More details on how to build  can be found in
Chapter [Chap:BuildingAMReX].

Running the Code
----------------

The example code can be run as follows,

::

      ./main3d.gnu.DEBUG.ex

The result may look like,

::

      Hello world from AMReX version 17.05-30-g5775aed933c4-dirty

The version string means the current commit 5775aed933c4 (note
that the first letter g in g577.. is not part of the hash)
is based on 17.05 with 30 additional commits and the work tree is dirty (i.e. there are uncommitted changes).

In the GNUmakefile there are compilation options for DEBUG
mode (less optimized code with more error checking), dimensionality,
compiler type, and flags to enable MPI and/or OpenMP parallelism.
If there are multiple instances of a parameter, the last instance
takes precedence.

Parallelization
---------------

Now let’s build with MPI by typing make USE\_MPI=TRUE (alternatively
you can set USE\_MPI=TRUE in the GNUmakefile). This
should make an executable named main3d.gnu.DEBUG.MPI.ex. Note
MPI in the file name. You can then run,

::

      mpiexec -n 4 ./main3d.gnu.DEBUG.MPI.ex

The result may look like,

::

      MPI initialized with 4 MPI processes
      Hello world from AMReX version 17.05-30-g5775aed933c4-dirty

If the compilation fails, you are referred to
Chapter [Chap:BuildingAMReX] on how to configure the build
system.

If you want to build with OpenMP, type make USE\_OMP=TRUE.
This should make an executable named main3d.gnu.DEBUG.OMP.ex. Note
OMP in the file name. Make sure the OMP\_NUM\_THREADS
environment variable is set on your system. You can then run,

::

      ./main3d.gnu.DEBUG.OMP.ex

The result may look like,

::

      OMP initialized with 4 OMP threads
      Hello world from AMReX version 17.06-287-g51875485fe51-dirty

Note that you can build with both USE\_MPI=TRUE and USE\_OMP=TRUE.
You can then run,

::

      mpiexec -n 2 ./main3d.gnu.DEBUG.MPI.OMP.ex

The result may look like,

::

      MPI initialized with 2 MPI processes
      OMP initialized with 4 OMP threads
      Hello world from AMReX version 17.06-287-g51875485fe51-dirty

Example: Heat Equation Solver
=============================

We now look at a more complicated example at
amrex/Tutorials/Basic/HeatEquation\_EX1\_C and show how simulation
results can be visualized. This example solves the heat equation,

.. math:: \frac{\partial\phi}{\partial t} = \nabla^2\phi

using forward Euler temporal integration on a periodic domain.
We could use a 5-point (in 2D) or 7-point (in 3D) stencil, but for demonstration
purposes we spatially discretize the PDE by first constructing fluxes on cell faces, e.g.,

.. math:: F_{i+\myhalf,j} = \frac{\phi_{i+1,j}-\phi_{i,j}}{\Delta x},

and then taking the divergence to update the cells,

.. math::

   \phi_{i,j}^{n+1} = \phi_{i,j}^n 
   + \frac{\Delta t}{\Delta x}\left(F_{i+\myhalf,j}-F_{i-\myhalf,j}\right)
   + \frac{\Delta t}{\Delta y}\left(F_{i,j+\myhalf}-F_{i,j-\myhalf}\right)

Don’t worry about the implementation details of the code.
You will be able to understand the code in this example after
Chapter [Chap:Basics].

Building and Running the Code
-----------------------------

To build a 2D executable, type make DIM=2. This will generate
an executable named main2d.gnu.ex. To run it, type,

::

      ./main2d.gnu.DEBUG.ex inputs_2d

Note that the command takes a file inputs\_2d. When the run
finishes, you will have a number of plotfiles, plt00000, plt01000, etc. The calculation solves the heat equation in 2D on a
:math:`256 \times 256` cells domain. It runs :math:`10,000` steps and makes a
plotfile every :math:`1,000` steps. These are runtime parameters that can
be adjusted in inputs\_2d.

Visualization
=============

There are several visualization tools that can be used for plotfiles. The standard tool used within the
-community is , a package developed and supported
by CCSE that is designed specifically for highly efficient visualization
of block-structured hierarchical AMR data.
Plotfiles can also be viewed using the , , and  packages.
Particle data can be viewed using .
Refer to Chapter [Chap:Visualization] for how to use each of these tools.
