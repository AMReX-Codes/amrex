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

The latest version of `yt` can be required for advanced options (e.g., rigid
injection for particles). To built `yt` directly from source, you can use

::

    pip install git+https://github.com/yt-project/yt.git


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

For a quick start-up, the most useful commands for post-processing can be found
in our Jupyter notebook
:download:`Visualization.ipynb<../../../Tools/Visualization.ipynb>`

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


Read back-transformed diagnotics
--------------------------------

When running a simulation in a boosted frame, WarpX has the capability to
back-transform the simulation results to the laboratory frame of reference, which
is often useful to study the physics. A set of function can be found in the
python file :download:`read_raw_data.py<../../../Tools/read_raw_data.py>`.

Alternatively, the main commands can be found in our example jupyter
notebook for postprocessing
:download:`Visualization.ipynb<../../../Tools/Visualization.ipynb>`.

The full back-transformed diagnostics of the entire domain is written in lab_frame_data/snapshots/ and the back-transformed diagnostics of the reduced domain is written to lab_frame_data/slices/
For instance : To plot the ``Ez`` field along the z-direction at the center of the 3D-domain of the full back-transformed diagnostics for the entire 3D domain :

::
    iteration = 0
    field = 'Ez'
    snapshot = './lab_frame_data/snapshots/' + 'snapshot' + str(iteration).zfill(5)
    header   = './lab_frame_data/snapshots/Header'
    allrd, info = read_raw_data.read_lab_snapshot(snapshot, header) # Read field data
    F = allrd[field]
    plt.plot(F[F.shape[0]//2,F.shape[1]//2-1,:])

Similarly, the back-transformed diagnostics on a reduced domain (1D line, 2D slice, 3D reduced diagnostic) can also be visualized using read_raw_data.py. For instance -- let us say that the user-input is an x-z slice at y=y_mid of the domain, then, to plot ``Ez`` of the x-z slice along the z-direction at the center of the slice :

::
    iteration = 0
    field = 'Ez'
    snapshot = './lab_frame_data/slices/' + 'slice' + str(iteration).zfill(5)
    header   = './lab_frame_data/slices/Header'
    allrd, info = read_raw_data.read_lab_snapshot(snapshot, header) # Read field data
    F_RD = allrd[field]
    plt.plot(F_RD[F_RD.shape[0]//2,0,:])


Note that, in the above snippet, we compare the 0th cell of the reduced diagnostic with F.shape[1]//2-1. For an x-z slice at y=y-mid of the domain, two cells are extracted at the center to ensure that the data format is HDF5 compliant. Let us consider that the domain consists of four cells in the y-dimension : [0,1,2,3], Then the 2D slice would contain the data that corresponds to [1,2]. That is the 0th cell of the reduced diagnostic corresponds to ny/2-1, (where, ny is the number of cells in the y-dimension).

If the back-transformed diagnostics are written in the HDF5 format (This can be done by compiling WarpX with USE_HDF5=TRUE), then the full domain snapshot and reduced domain diagnostics can be visualized using h5py :

::
    import matplotlib.pyplot as plt
    import h5py
    f1 = h5py.File('lab_frame_data/snapshots/snapshot00000', 'r')
    nx1 = f1['Ez'].shape[0]
    ny1 = f1['Ez'].shape[1]
    nz1 = f1['Ez'].shape[2]
    plt.plot(f1['Ez'][nx1//2,ny1//2-1,:])

    f2 = h5py.File('lab_frame_data/slices/slice00000', 'r')
    nx2 = f2['Ez'].shape[0]
    ny2 = f2['Ez'].shape[1]
    nz2 = f2['Ez'].shape[2]
    plt.figure()
    plt.plot(f2['Ez'][nx2//2,0,:])

The back-transformed particle data on the full and reduced diagnostic can be visualized as follows

::
    species='ions'
    iteration = 1

    snapshot = './lab_frame_data/snapshots/' + 'snapshot' + str(iteration).zfill(5)
    xbo = get_particle_field(snapshot, species, 'x') # Read particle data
    ybo = get_particle_field(snapshot, species, 'y')
    zbo = get_particle_field(snapshot, species, 'z')

    snapshot = './lab_frame_data/slices/' + 'slice' + str(iteration).zfill(5)
    xbo_slice = get_particle_field(snapshot, species, 'x') # Read particle data
    ybo_slice = get_particle_field(snapshot, species, 'y')
    zbo_slice = get_particle_field(snapshot, species, 'z')
    plt.figure()
    plt.plot(xbo, ybo, 'r.', markersize=1.)
    plt.plot(xbo_slice, ybo_slice, 'bx', markersize=1.)

Further information
-------------------

A lot more information can be obtained from the yt documentation, and the
corresponding notebook tutorials `here <http://yt-project.org/doc/>`__.
