Visualizing the simulation results
==================================

WarpX can write data either in `plotfile` format (AMReX's native format), or
in `openPMD format <https://www.openpmd.org/>`_ (a common data format for
Particle-In-Cell codes).

.. note::

    This is controlled by the parameters ``warpx.dump_plotfiles`` and
    ``warpx.dump_openpmd`` in the section :doc:`../running_cpp/parameters`.

This section describes some of the tools available to visualize the data:

.. toctree::
   :maxdepth: 1

   yt
   visit
   picviewer
   openpmdviewer
   advanced


In addition, WarpX also has In-Situ Visualization capabilities (i.e.
visualizing the data directly from the simulation, without dumping data
files to disk).

.. toctree::
    :maxdepth: 1

    sensei
    ascent
