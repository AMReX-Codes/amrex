#!/bin/bash
set -e # Exit with nonzero exit code if anything fails

# Build Doxygen
echo "Build the Doxygen documentation"
cd Doxygen
doxygen doxygen.conf &> doxygen.out
# Doxygen exits 0 even when it emits warnings or errors (e.g. LaTeX failures
# while rendering \f$...\f$ formulas), so the exit status is not a reliable
# gate; scan the log instead. Match both "warning:" (documentation problems)
# and "error:" (e.g. "error: Problems running latex...").
failures=$(grep -nE "warning:|error:" doxygen.out || true)
if [ -n "$failures" ]; then
    echo "Doxygen warnings/errors detected! Failing..."
    echo "$failures"
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
