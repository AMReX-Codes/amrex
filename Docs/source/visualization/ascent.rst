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
:code:`insitu.start` to the first time step where in situ analysis should be
performed.

Visualization/Analysis Pipeline Configuration
---------------------------------------------
Ascent uses the file :code:`ascent_actions.yaml` to configure analysis and
visualization pipelines. Ascent looks for the :code:`ascent_actions.yaml` file
in the current working directory.

For example, the following :code:`ascent_actions.yaml`
file extracts an isosurface of the field Ex for 15 levels and saves the
resulting images to :code:`levels_<nnnn>.png`. `Ascent Actions
<https://ascent.readthedocs.io/en/latest/Actions/index.html>`_ provides an
overview over all available analysis and visualization actions.

.. code-block:: json

    -
      action: "add_pipelines"
      pipelines:
        p1:
          f1:
            type: "contour"
            params:
               field: "Ex"
               levels: 15
    -
      action: "add_scenes"
      scenes:
        scene1:
          image_prefix: "levels_%04d"
          plots:
            plot1:
              type: "pseudocolor"
              pipeline: "p1"
              field: "Ex"

Here is another :code:`ascent_actions.yaml` example that renders isosurfaces
and particles:

.. code-block:: json

    -
      action: "add_pipelines"
      pipelines:
        p1:
          f1:
            type: "contour"
            params:
               field: "Bx"
               levels: 3
    -
      action: "add_scenes"
      scenes:
        scene1:
          plots:
            plot1:
              type: "pseudocolor"
              pipeline: "p1"
              field: "Bx"
            plot2:
              type: "pseudocolor"
              field: "particle_electrons_Bx"
              points:
                radius: 0.0000005
          renders:
            r1:
              camera:
                azimuth: 100
                elevation: 10
              image_prefix: "out_render_3d_%06d"


Finally, here is a more complex :code:`ascent_actions.yaml` example that
creates the same images as the prior example, but adds a trigger that
creates a Cinema Database at cycle 300:

.. code-block:: json

    -
      action: "add_triggers"
      triggers:
        t1:
          params:
            condition: "cycle() == 300"
            actions_file: "trigger.yaml"
    -
      action: "add_pipelines"
      pipelines:
        p1:
          f1:
            type: "contour"
            params:
               field: "jy"
               iso_values: [ 1000000000000.0, -1000000000000.0]
    -
      action: "add_scenes"
      scenes:
        scene1:
          plots:
            plot1:
              type: "pseudocolor"
              pipeline: "p1"
              field: "jy"
            plot2:
              type: "pseudocolor"
              field: "particle_electrons_w"
              points:
                radius: 0.0000002
          renders:
            r1:
              camera:
                azimuth: 100
                elevation: 10
              image_prefix: "out_render_jy_part_w_3d_%06d"


When the trigger condition is meet, `cycle() == 300`, the actions in
:code:`trigger.yaml` are also executed:

.. code-block:: json

    -
      action: "add_pipelines"
      pipelines:
        p1:
          f1:
            type: "contour"
            params:
               field: "jy"
               iso_values: [ 1000000000000.0, -1000000000000.0]
    -
      action: "add_scenes"
      scenes:
        scene1:
          plots:
            plot1:
              type: "pseudocolor"
              pipeline: "p1"
              field: "jy"
            plot2:
              type: "pseudocolor"
              field: "particle_electrons_w"
              points:
                radius: 0.0000001
          renders:
            r1:
              type: "cinema"
              phi: 10
              theta: 10
              db_name: "cinema_out"

You can view the Cinema Database result by opening
:code:`cinema_databases/cinema_out/index.html`.
