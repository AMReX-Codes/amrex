Running WarpX from PICMI
========================

See `PICMI information and source code <https://github.com/picmi-standard/picmi>`__.

Generate WarpX input from PICMI python script
---------------------------------------------

You can **install** the PICMI standard via, for example:

::

    pip install picmistandard

For other options, please visit link above.

To install the WarpX implementation of PICMI go to the WarpX directory and do:

::

    cd Python
    python setup.py install

Then you can **run**, for example, the :download:`Example PICMI input script for laser-plasma acceleration <../../../Examples/Physics_applications/laser_acceleration/PICMI_inputs_laser_acceleration.py>` with

::

    python PICMI_inputs_laser_acceleration.py

to create the WarpX input file that you can run with the WarpX executable.

Running WarpX with python binding form PICMI
--------------------------------------------

This option is not (yet) supported.
