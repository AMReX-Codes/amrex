In situ Visualization with SENSEI
=================================
SENSEI is a light weight framework for in situ data analysis. SENSEI's data
model and API provide uniform access to and run time selection of a diverse set
of visualization and analysis back ends including VisIt Libsim, ParaView
Catalyst, VTK-m, Ascent, ADIOS, Yt, and Python.

SENSEI uses an XML file to select and configure one or more back ends at run
time. Run time selection of the back end via XML means one user can access
Catalyst, another Libsim, yet another Python with no changes to the code.

Compiling with GNU Make
-----------------------
For codes making use of AMReX's build system add the following variable to the
code's main :code:`GNUmakefile`.

.. code-block:: bash

   USE_SENSEI_INSITU = TRUE

When set, AMReX's make files will query environment variables for the lists of
compiler and linker flags, include directories, and link libraries. These lists
can be quite elaborate when using more sophisticated back ends, and are best
set automatically using the :code:`sensei_config` command line tool that should
be installed with SENSEI. Prior to invoking make use the following command to
set these variables:

.. code-block:: bash

   source sensei_config

Typically, the :code:`sensei_config` tool is in the users PATH after loading
the desired SENSEI module. After configuring the build environment with
:code:`sensei_config`, proceed as usual.

.. code-block:: bash

   make -j4 -f GNUmakefile

ParmParse Configuration
-----------------------
Once an AMReX code has been compiled with SENSEI features enabled, it will need
to be enabled and configured at runtime. This is done using ParmParse input file.
The following 3 ParmParse parameters are used:

.. code-block:: python

   insitu.int = 2
   insitu.start = 0
   insitu.config = render_iso_catalyst_2d.xml

:code:`insitu.int` turns in situ processing on or off and controls how often
data is processed. :code:`insitu.start` controls when in situ processing
starts. :code:`insitu.config` points to the SENSEI XML file which selects and
configures the desired back end.

Obtaining SENSEI
-----------------
SENSEI is hosted on Kitware's Gitlab site at https://gitlab.kitware.com/sensei/sensei
It's best to checkout the latest release rather than working on the master branch.

To ease the burden of wrangling back end installs SENSEI provides two platforms
with all dependencies pre-installed, a VirtualBox VM, and a NERSC Cori
deployment. New users are encouraged to experiment with one of these.


SENSEI VM
~~~~~~~~~
The SENSEI VM comes with all of SENSEI's dependencies and the major back ends
such as VisIt and ParaView installed. The VM is the easiest way to test things
out. It also can be used to see how installs were done and the environment
configured.

NERSC Cori
~~~~~~~~~~
SENSEI is deployed at NERSC on Cori. The NERSC deployment includes the major
back ends such as ParaView Catalyst, VisIt Libsim, and Python.
