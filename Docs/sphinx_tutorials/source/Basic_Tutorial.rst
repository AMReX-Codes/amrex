.. role:: cpp(code)
   :language: c++

.. role:: fortran(code)
   :language: fortran

Tutorials/Basic
==========================

The tutorials in amrex/Tutorials/Basic demonstrate the most fundamental 
operations supported by AMReX.

For example, 

HelloWorld_C and HellowWorld_F demonstrate the GNU Make system -- with
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


