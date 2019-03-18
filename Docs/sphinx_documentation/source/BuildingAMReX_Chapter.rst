.. _Chap:BuildingAMReX:

Building AMReX
===================

In this chapter, we discuss AMReX's build systems. There are three ways to use
AMReX. Most AMReX developers use GNU Make. With this approach, there is no
installation step;  application codes adopt AMReX's build system and compile
AMReX while compiling their own codes. This will be discussed in more detail in
the section on :ref:`sec:build:make`.  The second approach is to build install
AMReX as a library (:ref:`sec:build:lib`); an application code then uses
its own build system and links to AMReX as an external library.  Finally, AMReX
can also be built with CMake, as detailed in the section on
:ref:`sec:build:cmake`.

Please note that we fully support AMReX for Linux systems in general and on the
DOE supercomputers (e.g. Cori, Summit) in particular.  Many of our users do build
and use AMReX on Macs but we do not have the resources to fully support Mac users.

.. toctree::
   :maxdepth: 2

   BuildingAMReX
