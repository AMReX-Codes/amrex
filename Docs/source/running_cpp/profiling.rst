Profiling the code
==================

Profiling with AMREX's built-in profiling tools
-----------------------------------------------
See `this page <https://amrex-codes.github.io/amrex/docs_html/Chapter12.html>`__ in the AMReX documentation.


Profiling the code with Intel advisor on NERSC
----------------------------------------------

Follow these steps:

- Instrument the code during compilation

  ::

     module swap craype-haswell craype-mic-knl
     make -j 16 COMP=intel USE_VTUNE=TRUE

   (where the first line is only needed for KNL)

- 

   
