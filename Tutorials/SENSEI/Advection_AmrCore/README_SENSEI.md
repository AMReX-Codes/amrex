# Using SingleVortex miniapp with SENSEI in situ #
SENSEI is a middleware that allows one to send data to various visualization and
analysis back ends through a uniform interface. It's data model and API enable
one to chose the desired visualization and analysis back end for a given task
with out limitting ones options, as the back ends can be inter-changed at run
time via a text based config file.

## Configuring the environment ##
First select the desired SENSEI install. Each install will support different set of
backends. This is necessary because not all of the back ends are compatible with
each other.

For instance to use SENSEI with ParaView Catalyst:
```
module load sensei/2.1.0-catalyst
```

To use SENSEI with VisIt Libsim:
```
module load sensei/2.1.0-libsim
```

## Compiling ##
SENSEI features in AMReX are conditionally compiled when the Make file variable
`USE_SENSEI_INSITU` is set. When this variable is set, the Make file will query
environment variables to determine the list of include directories and link
libraries needed to compile with SENSEI.

The `sensei_config` tool that is installed with SENSEI  will set the environment
variables that are used in the make files.

With a SENSEI module loaded, in a bash shell:
```
source sensei_config
```

Edit the file `Exec/SingleVortex/GNUmakefile`, set
```
USE_SENSEI_INSITU = TRUE
```

Finally, make the tutorial
```
make -j4
```

## Running ##
To use SENSEI in AMReX one needs to enable it via ParmParse input file.
Additionally one needs to provide a SENSEI XML configuration that selects
and configures the desired SENSEI backend.

Example XML configs are included in `Exec/SingleVortex/SENSEI`.

Edit the file `Exec/SingleVortex/inputs`

Running with ParaView Catalyst:
```
sensei.enabled = 1                          # turn SENSEI in situ on/off
sensei.config = SENSEI/render_catalyst.xml  # render simulation data with ParaView Catalyst
sensei.frequency = 1                        # number of level 0 steps between in situ processing
```

Running with VisIt Libsim:
```
sensei.enabled = 1                          # turn SENSEI in situ on/off
sensei.config = SENSEI/render_libsim.xml    # render simulation data with VisIt Libsim
sensei.frequency = 1                        # number of level 0 steps between in situ processing
```

Once the inputs files has been edited, run the execcutable as usual
```
mpiexec -np 4 ./main2d.gnu.MPI.ex inputs
```
