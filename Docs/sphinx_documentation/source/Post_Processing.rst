.. role:: cpp(code)
   :language: c++

.. _sec:amrvis:

WritePlotfileToASCII
====================

This basic routine reads in a single-level plotfile and writes the entire contents
to the standard output, one line at a time for each data value.
After reading in the plotfile to a :cpp:`MultiFab`, the program copies the data
into a separate :cpp:`MultiFab` with one large grid to make writing the data out
sequentially an easier task.

In ``amrex/Tools/Postprocessing/C_Src``, edit ``GNUMakefile`` to read
``EBASE = WritePlotfileToASCII`` and ``NEEDS_f90_SRC = FALSE`` and then ``make``
to generate an executable.  To run the executable, ``<executable> infile=<plotfilename>``.
You can modify the cpp file to write out on certain components, coordinates,
row/column formatting, etc.
