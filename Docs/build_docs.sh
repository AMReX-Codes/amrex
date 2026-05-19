#!/bin/bash
set -e # Exit with nonzero exit code if anything fails

# Build Doxygen
echo "Build the Doxygen documentation"
cd Doxygen
doxygen doxygen.conf &> doxygen.out
if grep -q "warning:" doxygen.out; then
    echo "Doxygen warnings detected! Failing..."
    cat doxygen.out
    exit 1
fi
cd ..

# Build sphinx
cd sphinx_documentation
echo "Build the Sphinx documentation"
make clean
make PYTHON="python3" latexpdf
cp build/latex/amrex.pdf source/
make SPHINXOPTS='-v -W --keep-going' PYTHON="python3" html
cd ..
