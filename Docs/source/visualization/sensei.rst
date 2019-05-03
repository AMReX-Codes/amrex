In situ Visualization with SENSEI
=================================
SENSEI is a light weight framework for in situ data analysis. SENSEI's data
model and API provide uniform access to and run time selection of a diverse set
of visualization and analysis back ends including VisIt Libsim, ParaView
Catalyst, VTK-m, Ascent, ADIOS, Yt, and Python.

SENSEI uses an XML file to select and configure one or more back ends at run
time. Run time selection of the back end via XML means one user can access
Catalyst, another Libsim, yet another Python with no changes to the code.

System Architecture
-------------------

.. _sensei_arch:
.. figure:: https://data.kitware.com/api/v1/item/5c06cd538d777f2179d4aaca/download

   SENSEI's in situ architecture enables use of a diverse of back ends which
   can be selected at run time via an XML configuration file

The three major architectural components in SENSEI are *data adaptors* which
present simulation data in SENSEI's data model, *analysis adaptors* which
present the back end data consumers to the simulation, and *bridge code* from
which the simulation manages adaptors and periodically pushes data through the
system. SENSEI comes equipped with a number of analysis adaptors enabling use
of popular analysis and visualization libraries such as VisIt Libsim, ParaView
Catalyst, Python, and ADIOS to name a few. AMReX contains SENSEI data adaptors
and bridge code making it easy to use in AMReX based simulation codes.

SENSEI provides a *configurable analysis adaptor* which uses an XML file to
select and configure one or more back ends at run time. Run time selection of
the back end via XML means one user can access Catalyst, another Libsim, yet
another Python with no changes to the code.  This is depicted in figure
:numref:`sensei_arch`. On the left side of the figure AMReX produces data, the
bridge code pushes the data through the configurable analysis adaptor to the
back end that was selected at run time.

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
The supported parameters are described in the following table.

+-------------------------+------------------------------------------------------+---------+
| parameter               | description                                          | default |
+=========================+======================================================+=========+
| :code:`insitu.int`      | turns in situ processing on or off and controls how  |    0    |
|                         | often data is processed.                             |         |
+-------------------------+------------------------------------------------------+---------+
| :code:`insitu.start`    | controls when in situ processing starts.             |    0    |
+-------------------------+------------------------------------------------------+---------+
| :code:`insitu.config`   | points to the SENSEI XML file which selects and      |         |
|                         | configures the desired back end.                     |         |
+-------------------------+------------------------------------------------------+---------+
| :code:`insitu.pin_mesh` | when 1 lower left corner of the mesh is pinned to    |    0    |
|                         | 0.,0.,0.                                             |         |
+-------------------------+------------------------------------------------------+---------+

A typical use case is to enabled SENSEI by setting :code:`insitu.int` to be
greater than 1, and :code:`insitu.config` to point SENSEI to an XML file that
selects and configures the desired back end.

.. code-block:: python

   insitu.int = 2
   insitu.config = render_iso_catalyst.xml

Back-end Selection and Configuration
------------------------------------
The back end is selected and configured at run time using the SENSEI XML file.
The XML sets parameters specific to SENSEI and to the chosen back end. Many of
the back ends have sophisticated configuration mechanisms which SENSEI makes
use of.  For example the following XML configuration was used on NERSC's Cori
with WarpX to render 10 iso surfaces, shown in figure :numref:`lpa_visit`, using
VisIt Libsim.

.. code-block:: xml

   <sensei>
     <analysis type="libsim" frequency="1" mode="batch"
       session="beam_j_pin.session"
       image-filename="beam_j_pin_%ts" image-width="1200" image-height="900"
       image-format="png" enabled="1"/>
   </sensei>

The *session* attribute names a session file that contains VisIt specific
runtime configuration. The session file is generated using VisIt GUI on a
representative dataset. Usually this data set is generated in a low resolution
run of the desired simulation.

.. _lpa_visit:
.. figure:: https://data.kitware.com/api/v1/item/5c06b4b18d777f2179d4784c/download

   Rendering of 10 3D iso-surfaces of j using VisIt libsim. The upper left
   quadrant has been clipped away to reveal innner structure.

The same run and visualization was repeated using ParaView Catalyst, shown in
figure :numref:`lpa_pv`, by providing the following XML configuration.

.. code-block:: xml

   <sensei>
     <analysis type="catalyst" pipeline="pythonscript"
       filename="beam_j.py" enabled="1" />
   </sensei>

Here the *filename* attribute is used to pass Catalyst a Catalyst specific
configuration that was generated using the ParaView GUI on a representative
dataset.

.. _lpa_pv:
.. figure:: https://data.kitware.com/api/v1/item/5c05b6388d777f2179d207ae/download

   Rendering of 10 3D iso-surfaces of j using ParaView Catalyst. The upper left
   quadrant has been clipped away to reveal innner structure.

The renderings in these runs were configured using a representative dataset
which was obtained by running the simulation for a few time steps at a lower
spatial resolution.  When using VisIt Libsim the following XML configures the
VTK writer to write the simulation data in VTK format. At the end of the run a
:code:`.visit` file that VisIt can open will be generated.

.. code-block:: xml

   <sensei>
     <analysis type="PosthocIO" mode="visit" writer="xml"
        ghost_array_name="avtGhostZones" output_dir="./"
        enabled="1">
     </analysis>
   </sensei>

When using ParaView Catalyst the following XML configures the VTK writer to
write the simulation data in VTK format. At the end of the run a :code:`.pvd`
file that ParaView can open will be generated.

.. code-block:: xml

   <sensei>
     <analysis type="PosthocIO" mode="paraview" writer="xml"
        ghost_array_name="vtkGhostType" output_dir="./"
        enabled="1">
     </analysis>
   </sensei>


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

The SENSEI VM can be downloaded here_.

.. _here: https://data.kitware.com/api/v1/file/5be656368d777f21799ee5a6/download

The SENSEI VM uses modules to manage the build and run environment. Load the
SENSEI modulefile for the back-end you wish to use. The following table
describes the available installs and which back-ends are supported in each.

+-------------------------------+-------------------------------------+
| modulefile                    | back-end(s)                         |
+===============================+=====================================+
| sensei/2.1.1-catalyst-shared  | ParaView Catalyst, ADIOS, Python    |
+-------------------------------+-------------------------------------+
| sensei/2.1.1-libsim-shared    | VisIt Libsim, ADIOS, Python         |
+-------------------------------+-------------------------------------+
| sensei/2.1.1-vtk-shared       | VTK-m, ADIOS, Python                |
+-------------------------------+-------------------------------------+

NERSC Cori
~~~~~~~~~~
SENSEI is deployed at NERSC on Cori. The NERSC deployment includes the major
back ends such as ADIOS, ParaView Catalyst, VisIt Libsim, and Python.

The SENSEI installs uses modules to manage the build and run environment. Load the
SENSEI modulefile for the back-end you wish to use. The following table
describes the available installs and which back-ends are supported in each.

+-------------------------------+-------------------------------------+
| modulefile                    | back-end(s)                         |
+===============================+=====================================+
| sensei/2.1.0-catalyst-shared  | ParaView Catalyst, ADIOS, Python    |
+-------------------------------+-------------------------------------+
| sensei/2.1.0-libsim-shared    | VisIt Libsim, ADIOS, Python         |
+-------------------------------+-------------------------------------+
| sensei/2.1.0-vtk-shared       | VTK-m, ADIOS, Python                |
+-------------------------------+-------------------------------------+


To access the SENSEI modulefiles on cori first add the SENSEI install to the search path:

.. code-block:: bash

    module use /usr/common/software/sensei/modulefiles


3D LPA Example
--------------
This section shows an example of using SENSEI and three different back ends on
a 3D LPA simulation. The instructions are specifically for NERSC cori, but also
work with the SENSEI VM. The primary difference between working through the examples
on cori or the VM are that different versions of software are installed.


Rendering with VisIt Libsim
~~~~~~~~~~~~~~~~~~~~~~~~~~~
First, log into cori and clone the git repo's.

.. code-block:: bash

   cd $SCRATCH
   mkdir warpx
   cd warpx/
   git clone --branch dev https://github.com/ECP-WarpX/WarpX.git WarpX-libsim
   git clone --branch development https://github.com/AMReX-Codes/amrex
   git clone --branch master https://bitbucket.org/berkeleylab/picsar.git
   cd WarpX-libsim
   vim GNUmakefile

Next, edit the makefile to turn the SENSEI features on.

.. code-block:: python

   USE_SENSEI_INSITU=TRUE

Then, load the SENSEI VisIt module, bring SENSEI's build requirements into the
environment, and compile WarpX.

.. code-block:: bash

   module use /usr/common/software/sensei/modulefiles/
   module load sensei/2.1.0-libsim-shared
   source sensei_config
   make -j8

Download the WarpX input deck, SENSEI XML configuration and and VisIt session
files. The inputs file configures WarpX, the xml file configures SENSEI, and
the session file configures VisIt. The inputs and xml files are written by
hand, while the session file is generated in VisIt gui on a representative data
set.

.. code-block:: bash

   wget https://data.kitware.com/api/v1/item/5c05d48e8d777f2179d22f20/download -O inputs.3d
   wget https://data.kitware.com/api/v1/item/5c05d4588d777f2179d22f16/download -O beam_j_pin.xml
   wget https://data.kitware.com/api/v1/item/5c05d4588d777f2179d22f0e/download -O beam_j_pin.session

To run the demo, submit an interactive job to the batch queue, and launch WarpX.

.. code-block:: bash

   salloc -C haswell -N 1 -t 00:30:00 -q debug
   ./Bin/main3d.gnu.TPROF.MPI.OMP.ex inputs.3d


Rendering with ParaView Catalyst
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
First, log into cori and clone the git repo's.

.. code-block:: bash

   cd $SCRATCH
   mkdir warpx
   cd warpx/
   git clone --branch dev https://github.com/ECP-WarpX/WarpX.git WarpX-catalyst
   git clone --branch development https://github.com/AMReX-Codes/amrex
   git clone --branch master https://bitbucket.org/berkeleylab/picsar.git
   cd WarpX-catalyst
   vim GNUmakefile

Next, edit the makefile to turn the SENSEI features on.

.. code-block:: python

   USE_SENSEI_INSITU=TRUE

Then, load the SENSEI ParaView module, bring SENSEI's build requirements into the
environment, and compile WarpX.

.. code-block:: bash

   module use /usr/common/software/sensei/modulefiles/
   module load sensei/2.1.0-catalyst-shared
   source sensei_config
   make -j8

Download the WarpX input deck, SENSEI XML configuration and and ParaView session
files. The inputs file configures WarpX, the xml file configures SENSEI, and
the session file configures ParaView. The inputs and xml files are written by
hand, while the session file is generated in ParaView gui on a representative data
set.

.. code-block:: bash

   wget https://data.kitware.com/api/v1/item/5c05b3fd8d777f2179d2067d/download -O inputs.3d
   wget https://data.kitware.com/api/v1/item/5c05b3fd8d777f2179d20675/download -O beam_j.xml
   wget https://data.kitware.com/api/v1/item/5c05b3fc8d777f2179d2066d/download -O beam_j.py

To run the demo, submit an interactive job to the batch queue, and launch WarpX.

.. code-block:: bash

   salloc -C haswell -N 1 -t 00:30:00 -q debug
   ./Bin/main3d.gnu.TPROF.MPI.OMP.ex inputs.3d

In situ Calculation with Python
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SENSEI's Python back-end loads a user provided script file containing callbacks
for :code:`Initialize`, :code:`Execute`, and :code:`Finalize` phases of the run.
During the execute phase the simulation pushes data through SENSEI.  SENSEI forwards
this data to the user provided Python function. SENSEI's MPI communicator is made
available to the user's function via a global variable :code:`comm`.

Here is a template for the user provided Python code.

.. code-block:: python

   # YOUR IMPORTS HERE

   # SET DEFAULTS OF GLOBAL VARIABLES THAT INFLUENCE RUNTIME BEHAVIOR HERE

   def Initialize():
     """ Initialization code """
     # YOUR CODE HERE
     return

   def Execute(dataAdaptor):
     """ Use sensei::DataAdaptor instance passed in
         dataAdaptor to access and process simulation data """
     # YOUR CODE HERE
     return

   def Finalize():
     """ Finalization code """
     # YOUR CODE HERE
     return

:code:`Initialize` and :code:`Finalize` are optional and will be called if
they are provided. :code:`Execute` is required. SENSEI's DataAdaptor API
is used to obtain data and metadata from the simulation. Data is through
VTK Object's. In WarpX the vtkOverlappingAMR VTK dataset is used.

The following script shows a simple integration of a scalar quantity
over the valid cells of the mesh. The result is saved in a CSV format.

.. code-block:: python

   import numpy as np, matplotlib.pyplot as plt
   from vtk.util.numpy_support import *
   from vtk import vtkDataObject
   import sys

   # default values of control parameters
   array = ''
   out_file = ''

   def Initialize():
     # rank zero writes the result
     if comm.Get_rank() == 0:
       fn = out_file if out_file else 'integrate_%s.csv'%(array)
       f = open(fn, 'w')
       f.write('# time, %s\n'%(array))
       f.close()
     return

   def Execute(adaptor):
     # get the mesh and arrays we need
     dobj = adaptor.GetMesh('mesh', False)
     adaptor.AddArray(dobj, 'mesh', vtkDataObject.CELL, array)
     adaptor.AddGhostCellsArray(dobj, 'mesh')
     time = adaptor.GetDataTime()

     # integrate over the local blocks
     varint = 0.
     it = dobj.NewIterator()
     while not it.IsDoneWithTraversal():
       # get the local data block and its props
       blk = it.GetCurrentDataObject()

       # get the array container
       atts = blk.GetCellData()

       # get the data array
       var =  vtk_to_numpy(atts.GetArray(array))

       # get ghost cell mask
       ghost = vtk_to_numpy(atts.GetArray('vtkGhostType'))
       ii = np.where(ghost == 0)[0]

       # integrate over valid cells
       varint = np.sum(var[ii])*np.prod(blk.GetSpacing())

       it.GoToNextItem()

     # reduce integral to rank 0
     varint = comm.reduce(varint, root=0, op=MPI.SUM)

     # rank zero writes the result
     if comm.Get_rank() == 0:
       fn = out_file if out_file else 'integrate_%s.csv'%(array)
       f = open(fn, 'a+')
       f.write('%s, %s\n'%(time, varint))
       f.close()
     return

The following XML configures SENSEI's Python back-end.

.. code-block:: xml

   <sensei>
     <analysis type="python" script_file="./integrate.py" enabled="1">
       <initialize_source>
   array='rho'
   out_file='rho.csv'
        </initialize_source>
     </analysis>
   </sensei>

The :code:`script_file` attribute sets the file path to load the user's Python
code from, and the :code:`initialize_source` element contains Python code that
controls runtime behavior specific to each user provided script.
