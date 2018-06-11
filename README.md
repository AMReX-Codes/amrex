# WarpX

Online documentation: https://ecp-warpx.github.io/index.html

## Overview

This repository is a project to build an **electromagnetic Particle-In-Cell (PIC) code** with mesh refinement, based on:

- The block-structure adaptive mesh refinement code [AMReX](https://bitbucket.org/berkeleylab/amrex)
- The highly-optimized PIC code PICSAR
- The Python interface of the code Warp

NB: [Regression tests](https://ccse.lbl.gov/pub/RegressionTesting/WarpX/) are run every night to check the validity of the code in this repository.

## Installation

### Downloading the code

You should clone the source codes of AMReX, PICSAR and WarpX into one single directory (e.g. `warpx_directory`):
```
mkdir warpx_directory
cd warpx_directory
git clone https://bitbucket.org/berkeleylab/warpx.git
git clone https://bitbucket.org/berkeleylab/picsar.git
git clone https://github.com/AMReX-Codes/amrex.git
```
You should then switch to the branch `development` of AMReX
```
cd amrex/
git checkout development
cd ..
```

### Compiling the code

`cd` into the the directory `warpx` and type
```
make -j 4
```
In this case, OpenMP is enabled by default. In order to disable OpenMP:
```
make -j 4 USE_OMP=FALSE
```
This will compile the code (using 4 procs for the compilation) and generate an executable file in the `Bin` directory.

In order to clean a previously compiled version:
```
make realclean
```

## Running the Langmuir tests

The folder `Example/Langmuir` contains code that allow the user
to run a Langmuir wave case with WarpX. The instructions below
explain how to do this.

### Running the test with Warpx

After compiling WarpX for (see the above instructions), copy the
compiled executable (its name starts with `main3d`) from the folder
`Bin/` to the folder
`Example/Langmuir/`. Then type
```
cd tests/Langmuir
mpirun -np 4 ./<executable_name> inputs
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
