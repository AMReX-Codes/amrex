.. role:: cpp(code)
   :language: c++

.. role:: fortran(code)
   :language: fortran

Tutorials/Basic
==========================

The tutorials in amrex/Tutorials/Basic demonstrate the most fundamental 
operations supported by AMReX.

**HelloWorld**
----------------

HelloWorld_C and HelloWorld_F demonstrate the GNU Make system -- with
a sample Make.package and GNUmakefile -- and the amrex::Initialize
and amrex::Finalize functions.

In addition, in HelloWorld_C, the amrex::Print() operation, 
which only prints from the I/O processor, is used to print out 
the AMReX version (as defined by amrex::Version()) being used. 

HelloWorld_F is a simple example of how to use the F_Interface routines,
which are Fortran wrappers for the underlying C++ data strutures and 
iterators.  Here, for example, rather than calling amrex::Print() in C++, we
test on whether amrex_parallel_ioprocessor() is true, and if so, invoke
the usual Fortran print call.

**main**
----------------

main_C and main_F introduce the following:

 1. By default, AMReX initializes MPI and uses MPI_COMM_WORLD as its communicator.
    However, applications could choose to initialize MPI themselves and pass in an
    existing communicator.

 2. By default, AMReX treats command line arguments as inputs parameters.  The expected
    format of argv is

        *executable inputs_file parm=value*

    Here, `executable` is the filename of the executable, `inputs_file` is the file containing
    runtime parameters used to build AMReX ParmParse database, and `parm=value` is an input
    parameter that will override its value in `inputs_file`.  Both `inputs_file` and
    `parm=value` are optional.  At most one `inputs_file` is allowed. Howeer, there can be
    multiple `parm=value` s.

    The parsing of the command line arguments is performed in amrex::Initialize.  Applications
    can choose to skip command line parsing.  Applications can also provide a function that
    adds parameters to AMReX ParmParse database.

**Build_with_libamrex**
-----------------------

This tutorial builds on the main_C example and demonstrates how to build the executable when we 
want to link local files with the pre-built amrex library (libamrex.a) that has been installed elsewhere.
We separate main.cpp from the main_C example into two separate files (main.cpp and 
test_parameters.cpp), replace MyAmr.H by MyParams.H and add a Fortran file my_func.f90.  
The GNUmakefile here assumes that you have already built the AMReX library; for instructions on how to do 
that see the  Building_libamrex_ chapter.

.. Building_libamrex: https://amrex-codes.github.io/amrex/docs_html/BuildingAMReX.html

**HeatEquation**
----------------

The HeatEquation examples solve a 2D or 3D (determined by how you set DIM in the GNUmakefile)
heat equation explicitly on a domain-decomposed mesh.  This example is described in detail in
the Basics_ chapter of the amrex Documentation

.. _Basics: https://amrex-codes.github.io/amrex/docs_html/Basics.html

