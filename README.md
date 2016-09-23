# WarpX

## Overview

This repository is a project to build an **electromagnetic Particle-In-Cell (PIC) code** with mesh refinement, based on:

- The block-structure adaptive mesh refinement code [BoxLib](http://boxlib-codes.github.io/)
- The highly-optimized PIC code PICSAR

## Installation

### Downloading the code

You should clone the source codes of Boxlib, PICSAR and WarpX into one single directory (e.g. `warpx_directory`):
```
mkdir warpx_directory
cd warpx_directory
git clone https://github.com/BoxLib-Codes/BoxLib.git
git clone https://bitbucket.org/berkeleylab/picsar.git
git clone https://bitbucket.org/remilehe/warpx.git
```
You should then switch to the branch `development` of BoxLib, and to the branch `picsar_mpi_com_opt` of picsar:
```
cd BoxLib/
git checkout development
cd ../picsar
git checkout picsar_mpi_com_opt
cd ..
```

**For Mac OSX users:**   
There is an additional step: open the file `BoxLib/Tools/C_mk/Make.mpi` and add at the **end of the file** 
```
LIBRARY_LOCATIONS += <your_MPI_library>
INCLUDE_LOCATIONS += <your_MPI_include>
```
where `<your_MPI_library>` and `<your_MPI_include>` should be replaced by the corresponding location on your system. For instance, if you installed openmpi and gcc4.8 with Macports, `<your_MPI_library>` should be `/opt/local/lib/openmpi-gcc48` and `<your_MPI_include>` should be `/opt/local/include/openmpi-gcc48`.

In addition, you may need to compile the PICSAR files with `mpif90`.


### Compiling the code

`cd` into the the directory `warpx` and type
```
make -j 4 USE_MPI=TRUE BOXLIB_USE_MPI_WRAPPERS=TRUE
```

In order to clean a previously compiled version:
```
make realclean
```

## Organization of the files

Several files in the directory `warpx` are key to the project:

- The file `main.cpp` contains the main program which is executed.
- The file `single_level.cpp` contains the function `single_level` which does most of the physical work.
- The file `Particles.H` was copied from BoxLib and modified to include PICSAR calls. The methods in `Particle.H` are used in `single_level`.
