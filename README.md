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
make -j
```

In order to clean a previously compiled version:
```
make realclean
```

## Running the tests

The folder `tests/Langmuir` contains code that allow the user
to run a Langmuir wave case with either WarpX or PICSAR, and to
compare the results. The instructions below explain how to do this.

### Running the test with Warpx

After compiling WarpX (see the above instructions), copy the
compiled executable (its name starts with `main3d`) to the folder
`tests/Langmuir/`. Then type
```
cd tests/Langmuir
./<executable_name> input_warpx
```
where `<executable_name>` should be replaced by the name of the
compiled executable.

The code produces a set of folders with names starting with `plt`.

### Running the test with PICSAR

`cd` into the folder `picsar` (`cd ../../../picsar` if you are
currently in the folder `test/Langmuir`), and type the following set
of commands:
```
make clean
make
cp fortran_bin/picsar ../warpx/tests/Langmuir/
cd ../warpx/tests/Langmuir/
./picsar
```

## Using WarpX

### Configuration of the input file

current_deposition_algo: algorithm for the current deposition
 - 3: Scalar classical current deposition
 - 2: Optimized classical current
 - 1: Esirkepov non optimized
 - 0: Esirkepov optimized

charge_deposition_algo:
 - 0: Optimized version
 - 1: Scalar version

field_gathering_algo:
 - 0: Optmized subroutines
 - 1: Scalar subroutines
 - 2: General order non-optimized version

particle_pusher_algo: algorithm for the particle pusher
 - 0: pusher of Boris
 - 1: pusher of J. L. Vay

### Visualizing and comparing the results

The results are compared using Python, inside a Jupyter notebook. In
order to be able to visualize the results, you need to first install
the proper visualization software:

- If you are using the Anaconda distribution of Python, type:
```
conda install jupyter numpy yt matplotlib
```

- Otherwise, type:
```
pip install jupyter numpy yt matplotlib
```

Then, within the folder `tests/Langmuir`, type:
```
jupyter notebook Visualization.ipynb 
```
and follow the instructions that will pop up in your browser.
