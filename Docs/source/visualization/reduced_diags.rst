Reduced diagnostics
===================

WarpX has optional reduced diagnostics, that typically return one value (e.g., particle energy) per timestep. A Python function to read them is available from module ``read_raw_data`` in ``WarpX/Tools/``:

.. code-block:: python

   from read_raw_data import read_reduced_diags
   filename = 'EF.txt'
   metadata, data = read_raw_data.read_reduced_diags( filename )
   # list available diagnostics
   data.keys()
   # Print total field energy on level 0
   data['total_lev0']
   # Print units for the total field energy on level 0
   metadata['units']['total_lev0']
