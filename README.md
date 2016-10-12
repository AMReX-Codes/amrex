# WarpX

## Overview

This repository is a project to build an **electromagnetic Particle-In-Cell (PIC) code** with mesh refinement, based on:

- The block-structure adaptive mesh refinement code [BoxLib](http://boxlib-codes.github.io/)
- The highly-optimized PIC code PICSAR
- The Python interface of the code Warp

## Installation

### Downloading the code

You should clone the source codes of Boxlib, PICSAR and WarpX into one single directory (e.g. `warpx_directory`):
```
mkdir warpx_directory
cd warpx_directory
git clone https://github.com/BoxLib-Codes/BoxLib.git
git clone https://bitbucket.org/berkeleylab/picsar.git
git clone https://bitbucket.org/berkeleylab/warpx.git
```
You should then switch to the branch `development` of BoxLib
```
cd BoxLib/
git checkout development
cd ..
```

**For Mac OSX users:**   
There is an additional step: create a file
`BoxLib/Tools/C_mk/Make.local` and include the following text
```
CXX = mpicxx
CC = mpicc
FC = mpif90
fc = mpif90
F90 = mpif90

LIBRARY_LOCATIONS += <your_MPI_library>
INCLUDE_LOCATIONS += <your_MPI_include>
BL_MPI_LIBS = <your_compilation_flags>
```

where `<your_MPI_library>` and `<your_MPI_include>` should be replaced
by the corresponding location on your system. For instance, if you
installed openmpi and gcc4.8 with Macports, `<your_MPI_library>`
should be `/opt/local/lib/openmpi-gcc48` and `<your_MPI_include>`
should be `/opt/local/include/openmpi-gcc48`.

In addition, `<your_compilation_flags>` should be replaced by the
compilation flags that you see when typing `mpif90 -show`. For openmpi
installed with Macports, this is `-lmpi_usempi -lmpi_mpifh -lmpi`.

### Compiling the code

`cd` into the the directory `warpx` and type
```
make -j 4 USE_MPI=TRUE BOXLIB_USE_MPI_WRAPPERS=TRUE
```

In order to clean a previously compiled version:
```
make realclean
```

