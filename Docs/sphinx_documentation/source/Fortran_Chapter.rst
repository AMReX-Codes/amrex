.. role:: cpp(code)
   :language: c++

.. _Chap:Fortran:

Fortran Interface
=================


The core of AMReX is written in C++. For Fortran users who want to write all of
their programs in Fortran, AMReX provides Fortran interfaces around most of
functionalities except for the :cpp:`AmrLevel` class (see the chapter on
:ref:`Chap:AmrLevel`) and particles (see the chapter on :ref:`Chap:Particles`).
We should not confuse the Fortran interface in this chapter with the Fortran
kernel functions called inside :cpp:`MFIter` loops in codes (see the section
on :ref:`sec:basics:fortran`). For the latter, Fortran is used in some sense as
a domain-specific language with native multi-dimensional arrays, whereas here
Fortran is used to drive the whole application code. In order to better
understand AMReX, Fortran interface users should read the rest of the documentation
except for the Chapters on :ref:`Chap:AmrLevel` & :ref:`Chap:Particles`.


.. toctree::
   :maxdepth: 1

   Fortran
