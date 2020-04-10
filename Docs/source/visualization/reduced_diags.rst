Reduced diagnostics
===================

WarpX has optional reduced diagnostics, that typically return one value (e.g., particle energy) per timestep.

A simple and quick way to read the data using python is

.. code-block:: python

    data = numpy.genfromtxt("filename.txt")

where ``data`` is a two dimensional array, ``data[i][j]`` gives the data in the ith row and the jth column.

A Python function to read the data is available from module ``read_raw_data`` in ``WarpX/Tools/PostProcessing/``:

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

In addition, for reduced diagnostic type ``ParticleHistogram``,
another Python function is available:

.. code-block:: python

    from read_raw_data import read_reduced_diags_histogram
    filename = 'velocity_distribution.txt'
    metadata_dict, data_dict, bin_value, bin_data = read_reduced_diags_histogram( filename )
    # 1-D array of the ith bin value
    bin_value[i]
    # 2-D array of the jth bin data at the ith time
    bin_data[i][j]

