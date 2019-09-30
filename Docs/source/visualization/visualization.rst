Visualizing the simulation results
==================================

WarpX can write data either in `plotfile` format (AMReX's native format), or
in `openPMD format <https://www.openpmd.org/>`_ (a common data format for
Particle-In-Cell codes).

.. note::

    This is controlled by the parameters ``warpx.dump_plotfiles`` and
    ``warpx.dump_openpmd`` & ``warpx.openpmd_backend`` in the section
    :doc:`../running_cpp/parameters`.

This section describes some of the tools available to visualize the data:

.. toctree::
   :maxdepth: 1

   yt
   visit
   picviewer
   openpmdviewer
   advanced
   plot_parallel


In addition, WarpX also has In-Situ Visualization capabilities (i.e.
visualizing the data directly from the simulation, without dumping data
files to disk).

.. toctree::
    :maxdepth: 1

    sensei
    ascent

If you like the 3D rendering of laser wakefield acceleration
on the WarpX documentation frontpage (which is
also the avatar of the ECP-WarpX organization), you can find the serial
analysis script :download:`video_yt.py<../../../Tools/video_yt.py>` as well
as a parallel analysis script
:download:`video_yt.py<../../../Tools/yt3d_mpi.py>` used to make a similar
rendering for a beam-driven wakefield simulation, running parallel.
