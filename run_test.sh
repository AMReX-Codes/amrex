#!/bin/bash

# Create test directory
rm -rf test_dir
mkdir test_dir
cd test_dir

# Copy WarpX into current test directory
mkdir warpx
cp -r ../* warpx

# Clone PICSAR and AMReX
git clone https://github.com/AMReX-Codes/amrex.git
cd amrex ; git checkout development ; cd ..
git clone https://bitbucket.org/berkeleylab/picsar.git

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
