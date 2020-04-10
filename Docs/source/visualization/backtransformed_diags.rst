Visualizing back-transformed diagnostics
========================================

When running a simulation in a boosted frame, WarpX has the capability to
back-transform the simulation results to the laboratory frame of reference, which
is often useful to study the physics. A set of functions can be found in the
python file :download:`read_raw_data.py<../../../Tools/PostProcessing/read_raw_data.py>`. The main commands can be found in our example jupyter notebook for postprocessing :download:`Visualization.ipynb<../../../Tools/PostProcessing/Visualization.ipynb>`.

The full back-transformed diagnostics of the entire domain is written in ``lab_frame_data/snapshots/`` and the back-transformed diagnostics of the reduced domain is written to ``lab_frame_data/slices/``
For instance: To plot the ``Ez`` field along the z-direction at the center of the 3D-domain of the full back-transformed diagnostics for the entire 3D domain:

::

    import read_raw_data
    import matplotlib.pyplot as plt

    iteration = 0
    field = 'Ez'
    snapshot = './lab_frame_data/snapshots/' + 'snapshot' + str(iteration).zfill(5)
    header   = './lab_frame_data/snapshots/Header'
    allrd, info = read_raw_data.read_lab_snapshot(snapshot, header) # Read field data
    F = allrd[field]
    plt.plot(F[F.shape[0]//2,F.shape[1]//2-1,:])

Similarly, the back-transformed diagnostics on a reduced domain (1D line, 2D slice, 3D reduced diagnostic) can also be visualized using read_raw_data.py. For instance -- let us say that the user-input is an "x-z" slice (at the center of the domain in the "y-direction"), then, to plot ``Ez`` on this x-z slice:

::

    iteration = 0
    field = 'Ez'
    snapshot = './lab_frame_data/slices/' + 'slice' + str(iteration).zfill(5)
    header   = './lab_frame_data/slices/Header'
    allrd, info = read_raw_data.read_lab_snapshot(snapshot, header) # Read field data
    F_RD = allrd[field]
    plt.plot(F_RD[F_RD.shape[0]//2,0,:])


Note that, in the above snippet, we compare the 0th cell of the reduced diagnostic with ``F.shape[1]//2-1``. For an x-z slice at y=y-mid of the domain, two cells are extracted at the center to ensure that the data format is HDF5 compliant. Let us consider that the domain consists of four cells in the y-dimension: [0,1,2,3], Then the 2D slice would contain the data that corresponds to [1,2]. That is the 0th cell of the reduced diagnostic corresponds to ``ny/2-1``, (where, ny is the number of cells in the y-dimension).

If the back-transformed diagnostics are written in the HDF5 format (This can be done by compiling WarpX with USE_HDF5=TRUE), then the full domain snapshot and reduced domain diagnostics can be visualized using h5py:

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
