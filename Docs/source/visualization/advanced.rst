Advanced yt visualization, for developers (for plotfiles)
=========================================================

This sections contains yt commands for advanced users. The Particle-In-Cell methods uses a
staggered grid (see :doc:`../theory/picsar_theory`), so that the x, y, and z components of the
electric and magnetic fields are all defined at different locations in space. Regular output
(see the :doc:`yt` page, or the notebook at ``WarpX/Tools/Visualization.ipynb`` for an example)
returns cell-centered data for convenience, which involves an additional operation. It is sometimes
useful to access the raw data directly. Furthermore,
the WarpX implementation for mesh refinement contains a number of grids for each level (coarse,
fine and auxilary, see
:doc:`../theory/warpx_theory` for more details), and it is sometimes useful to access each of
them (regular output return the auxiliary grid only). This page provides information to read
raw data of all grids.

Dump additional data
--------------------

In order to dump additional data in WarpX (mostly for debugging purpose), run the simulation
with parameters

::

    warpx.plot_raw_fields = 1
    warpx.plot_finepatch = 1
    warpx.plot_crsepatch = 1
    warpx.plot_dive = 1
    warpx.plot_rho = 1

see :doc:`../running_cpp/parameters` for more information on these parameters.

Read raw data
-------------

Meta-data
relevant to this topic (number and locations of grids in the simulation) are accessed to
with

::

    import yt
    # get yt dataset
    ds = yt.load( './plotfiles/plt00004' )
    # Index of data in the plotfile
    ds_index = ds.index
    # Print the number of grids in the simulation
    ds_index.grids.shape
    # Left and right physical boundary of each grid
    ds_index.grid_left_edge
    ds_index.grid_right_edge
    # List available fields
    ds.field_list

When ``warpx.plot_raw_fields=1`` and ``warpx.plot_finepatch=1``, here are some useful
commands to access properties of a grid and the Ex field on the fine patch:

::

    # store grid number 2 into my_grid
    my_grid = ds.index.grids[2]
    # Get left and right edges of my_grid
    my_grid.LeftEdge
    my_grid.RightEdge
    # Get Level of my_grid
    my_grid.Level
    # left edge of the grid, in number of points
    my_grid.start_index

Return the ``Ex`` field on the fine patch of grid ``my_grid``:

::

    my_field = my_grid['raw', 'Ex_fp'].squeeze().v

For a 2D plotfile, ``my_field`` has shape ``(nx,nz,2)``. The last component stands for the
two values on the edges of each cell for the electric field, due to field staggering. Numpy
function ``squeeze`` removes empty components. While ``yt`` arrays are unit-aware, it is
sometimes useful to extract the data into unitless numpy arrays. This is achieved with ``.v``.
In the case of ``Ex_fp``, the staggering is on direction ``x``, so that
``my_field[:,:-1,1] == my_field[:,1:,0]``.

All combinations of the fields (``E`` or ``B``), the component (``x``, ``y`` or ``z``) and the
grid (``_fp`` for fine, ``_cp`` for coarse and ``_aux`` for auxiliary) can be accessed in this
way, i.e., ``my_grid['raw', 'Ey_aux']`` or ``my_grid['raw', 'Bz_cp']`` are valid queries.
