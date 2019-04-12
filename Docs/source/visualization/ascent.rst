In situ Visualization with Ascent
=================================
Ascent is a system designed to meet the in-situ visualization and analysis
needs of simulation code teams running multi-physics calculations on many-core
HPC architectures. It provides rendering runtimes that can leverage many-core
CPUs and GPUs to render images of simulation meshes.

Compiling with GNU Make
-----------------------
After building and installing Ascent according to the instructions at
`Building Ascent <https://ascent.readthedocs.io/en/latest/BuildingAscent.html>`_,
you can enable it in WarpX by changing the line

.. code-block:: bash

   USE_ASCENT_INSITU = FALSE

in GNUmakefile to 

.. code-block:: bash

   USE_ASCENT_INSITU = TRUE

Furthermore, you must ensure that either the :code:`ASCENT_HOME` shell
environment variable contains the directory where Ascent is installed
or you must specify this location when invoking make, i.e.,

.. code-block:: bash

   make -j 8 ASCENT_HOME = /path/to/ascent/install

ParmParse Configuration
-----------------------
Once an AMReX code has been compiled with Ascent enabled, it will need
to be enabled and configured at runtime. This is done using ParmParse input file.
The supported parameters are described in the following table.

+-------------------------+------------------------------------------------------+---------+
| parameter               | description                                          | default |
+=========================+======================================================+=========+
| :code:`insitu.int`      | turns in situ processing on or off and controls how  |    0    |
|                         | often data is processed.                             |         |
+-------------------------+------------------------------------------------------+---------+
| :code:`insitu.start`    | controls when in situ processing starts.             |    0    |
+-------------------------+------------------------------------------------------+---------+

A typical use case is setting :code:`insitu.int` to a value of one or greater and
:code:`insitu.start` to the first time step where in situ analyswhere in situ analysis should be
performed.

Visualization/Analysis Pipeline Configuration
---------------------------------------------
Ascent uses the file :code:`ascent_actions.json` to configure analysis and
visualization pipelines. For example, the following :code:`ascent_actions.json`
file extracts an isosurface of the field Ex for 15 levels and saves the
resulting images to :code:`levels_<nnnn>.png`. `Ascent Actions 
<https://ascent.readthedocs.io/en/latest/Actions/index.html>`_ provides an
overview over all available analysis and visualization actions.

.. code-block:: json

  [
    {
      "action": "add_pipelines",
      "pipelines": 
      {
        "p1": 
        {
          "f1": 
          {
            "type" : "contour",
            "params" :
            {
              "field" : "Ex",
              "levels": 15
            }
          }
        }
      }
    },
    {
      "action": "add_scenes",
      "scenes":
      {
        "s1": 
        {
          "image_prefix": "levels_%04d",
          "plots": 
          {
          "p1": 
            {
              "type": "pseudocolor",
              "pipeline": "p1",
              "field": "Ex"
            }
          }
        }
      }
    },
  
    {
      "action": "execute"
    },
  
    {
      "action": "reset"
    }
  ]

