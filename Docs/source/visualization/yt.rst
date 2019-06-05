Visualization with yt (for plotfiles)
=====================================

`yt <http://yt-project.org/>`__ is a Python package that can help in analyzing
and visualizing WarpX data (among other data formats). It is convenient
to use yt within a `Jupyter notebook <http://jupyter.org/>`__.

Installation
------------

From the terminal:

::

    pip install yt jupyter

or with the `Anaconda distribution <https://anaconda.org/>`__ of python (recommended):

::

    conda install -c conda-forge yt

Visualizing the data
--------------------

Once data ("plotfiles") has been created by the simulation, open a Jupyter notebook from
the terminal:

::

    jupyter notebook

Then use the following commands in the first cell of the notebook to import yt
and load the first plot file:

::

    import yt
    ds = yt.load('./diags/plotfiles/plt00000/')

The list of field data and particle data stored can be seen with:

::

    ds.field_list

Field data
~~~~~~~~~~

Field data can be visualized using ``yt.SlicePlot`` (see the docstring of
this function `here <http://yt-project.org/doc/reference/api/yt.visualization.plot_window.html#yt.visualization.plot_window.SlicePlot>`__)

For instance, in order to plot the field ``Ex`` in a slice orthogonal to ``y`` (``1``):

::

    yt.SlicePlot( ds, 1, 'Ex', origin='native' )

.. note::

    `yt.SlicePlot` creates a 2D plot with the same aspect ratio as the physical
    size of the simulation box. Sometimes this can lead to very elongated plots
    that are difficult to read. You can modify the aspect ratio with the
    `aspect` argument ; for instance:

    ::

        yt.SlicePlot( ds, 1, 'Ex', aspect=1./10 )


Alternatively, the data can be obtained as a `numpy <http://www.numpy.org/>`__ array.

For instance, in order to obtain the field `jz` (on level 0) as a numpy array:

::

    ad0 = ds.covering_grid(level=0, left_edge=ds.domain_left_edge, dims=ds.domain_dimensions)
    jz_array = ad0['jz'].to_ndarray()


Particle data
~~~~~~~~~~~~~

Particle data can be visualized using ``yt.ParticlePhasePlot`` (see the docstring
`here <http://yt-project.org/doc/reference/api/yt.visualization.particle_plots.html?highlight=particlephaseplot#yt.visualization.particle_plots.ParticlePhasePlot>`__).

For instance, in order to plot the particles' ``x`` and ``y`` positions:

::

    yt.ParticlePhasePlot( ds.all_data(), 'particle_position_x', 'particle_position_y', 'particle_weight')

Alternatively, the data can be obtained as a `numpy <http://www.numpy.org/>`__ array.

For instance, in order to obtain the array of position `x` as a numpy array:

::

    ad = ds.all_data()
    x = ad['particle_position_x'].to_ndarray()


Further information
-------------------

A lot more information can be obtained from the yt documentation, and the
corresponding notebook tutorials `here <http://yt-project.org/doc/>`__.

A notebook with useful yt commands can be found at ``WarpX/Tools/Visualization.ipynb``.
It is a good start to explore simulation results.
