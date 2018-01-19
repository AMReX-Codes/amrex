.. _BuildingAMReX:

Building AMReX
===================

In this chapter, we discuss AMReX's build systems. There are three
ways to use AMReX. The approach used by AMReX developers uses GNU
Make. There is no installation step in this approach. Application
codes adopt AMReX's build system and compile AMReX while compiling
their own codes. This will be discussed in more detail in the section 
on :ref:`sec:build:make`. The second approach is to build AMReX into
a library and installing it (cf. :ref:`sec:build:lib`). Then an
application code uses its own build system and links AMReX as an
external library. AMReX can also be built with Cmake, as detailed in 
the section on :ref:`sec:build:cmake`.



.. toctree::
   :maxdepth: 2

   BuildingAMReX
