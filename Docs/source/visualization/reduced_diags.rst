Reduced diagnostics
===================

WarpX has optional reduced diagnostics, that typically return one value (e.g., particle energy) per timestep.

A simple and quick way to read the data using python is

.. code-block:: python

    data = numpy.genfromtxt("filename.txt")

where ``data`` is a two dimensional array, ``data[i][j]`` gives the data in the ith row and the jth column.

In addition, a Python function to read the data is available from module ``read_raw_data`` in ``WarpX/Tools/``:

.. code-block:: python

   from read_raw_data import read_reduced_diags
   filename = 'EF.txt'
   metadata, data = read_reduced_diags( filename )
   # list available diagnostics
   data.keys()
   # Print total field energy on level 0
   data['total_lev0']
   # Print units for the total field energy on level 0
   metadata['units']['total_lev0']
