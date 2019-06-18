Visualization with openPMD-viewer (for openPMD data)
====================================================

openPMD-viewer is an open-source Python package to access openPMD data.

It allows to:
- Quickly browse through the data, with a GUI-type interface in the Jupyter notebook
- Have access to the data numpy array, for more detailed analysis

Installation
------------

openPMD-viewer can be installed via ``conda`` or ``pip``:

::

    conda install -c rlehe openpmd_viewer

::

    pip install openPMD-viewer

Usage
-----

openPMD-viewer can be used either in simple Python scripts, or in a Jupyter
notebook. In both cases, you can import openPMD-viewer, and load the data
with the following commands:

::

    from opmd_viewer import OpenPMDTimeSeries
    ts = OpenPMDTimeSeries('./diags/hdf5')

.. note::

    If you are using the Jupyter notebook, then you can start a pre-filled
    notebook, which already contains the above lines, by typing in a terminal:

    ::

        openPMD_notebook

When using the Jupyter notebook, you can quickly browse through the data
by using the command:

::

    ts.slider()

You can also access the particle and field data as numpy arrays with the
methods ``ts.get_field`` and ``ts.get_particle``. See the openPMD-viewer
tutorials `here <https://github.com/openPMD/openPMD-viewer/tree/master/tutorials>`_ for more info.
