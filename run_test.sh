#!/bin/bash
# This script runs some of WarpX's standard regression tests, but
# without comparing the output to previously run simulations.
# This checks that:
# - The code compiles and runs without error
# - For some of the tests, a Python script checks that the results are
# physically correct.

# The tests can be influenced by environment variables:
# Use `export WARPX_TEST_DIM=3` or `export WARPX_TEST_DIM=2` in order to
# select only the tests that correspond to this dimension
# Use `export WARPX_TEST_ARCH=CPU` or `export WARPX_TEST_ARCH=GPU` in order
# to run the tests on CPU or GPU respectively.

# Create test directory
rm -rf test_dir
mkdir test_dir
cd test_dir

# Copy WarpX into current test directory
mkdir warpx
cp -r ../* warpx

# Clone PICSAR and AMReX
git clone --branch development https://github.com/AMReX-Codes/amrex.git
git clone --branch master https://bitbucket.org/berkeleylab/picsar.git

# Clone the AMReX regression test utility
git clone https://github.com/RemiLehe/regression_testing.git

# Prepare regression tests
mkdir -p rt-WarpX/WarpX-benchmarks
cd warpx/Regression
python prepare_file_travis.py
cp travis-tests.ini ../../rt-WarpX

# Run the tests
cd ../../regression_testing/
python regtest.py ../rt-WarpX/travis-tests.ini --no_update all --source_git_hash=$WARPX_TEST_COMMIT
