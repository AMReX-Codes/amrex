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

# Remove contents and link to a previous test directory (intentionally two arguments)
rm -rf test_dir/* test_dir
# Create a temporary test directory
tmp_dir=$(mktemp --help >/dev/null 2>&1 && mktemp -d -t ci-XXXXXXXXXX || mktemp -d "${TMPDIR:-/tmp}"/ci-XXXXXXXXXX)
if [ $? -ne 0 ]; then
    echo "Cannot create temporary directory"
    exit 1
fi

# Copy WarpX into current test directory
mkdir ${tmp_dir}/warpx
cp -r ./* ${tmp_dir}/warpx

# Link the test directory
ln -s ${tmp_dir} test_dir

# Switch to the test directory
cd test_dir

# Clone PICSAR and AMReX
git clone --branch development https://github.com/AMReX-Codes/amrex.git
# Use QED brach for QED tests
if [ ${HAS_QED} = "TRUE" ]; then
    git clone --branch QED https://bitbucket.org/berkeleylab/picsar.git
else
    git clone --branch master https://bitbucket.org/berkeleylab/picsar.git
fi

# Clone the AMReX regression test utility
git clone https://github.com/ECP-WarpX/regression_testing.git

# Prepare regression tests
mkdir -p rt-WarpX/WarpX-benchmarks
cd warpx/Regression
python prepare_file_travis.py
cp travis-tests.ini ../../rt-WarpX

# Run the tests
cd ../../regression_testing/
python regtest.py ../rt-WarpX/travis-tests.ini --no_update all --source_git_hash=${WARPX_TEST_COMMIT}
