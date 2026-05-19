#!/bin/bash
set -e # Exit with nonzero exit code if anything fails

# Build Doxygen
echo "Build the Doxygen documentation"
cd Docs/Doxygen
doxygen doxygen.conf &> doxygen.out
if grep -q "warning:" doxygen.out; then
    echo "Doxygen warnings detected! Failing..."
    cat doxygen.out
    exit 1
fi
cd ../..

# remove LaTeX temp files left by doxygen formula rendering
cd Docs/Doxygen/html
rm -f *.repository *.aux *.dvi *.log *.tex *.ps *.epsi
cd ../../..

# copy doxygen to target location
mkdir build
cd build
mkdir docs_html docs_xml
mkdir -p docs_html/doxygen
cp -rp ../Docs/Doxygen/html/* docs_html/doxygen/
mkdir -p docs_xml/doxygen
# add tagfile to allow other docs to interlink with amrex.
# but we don't need to copy other xml files.
cp ../Docs/Doxygen/amrex-doxygen-web.tag.xml docs_xml/doxygen/.
cd ..

# Build sphinx
cd Docs/sphinx_documentation
echo "Build the Sphinx documentation for Amrex."
# Stabilize PDF metadata (CreationDate, ModDate, /ID) so identical
# content produces byte-identical output across repeated builds.
export SOURCE_DATE_EPOCH=$(date -d "$(date +%Y)-04-01 00:00:00 UTC" +%s) # yyyy-04-01 UTC
make PYTHON="python3" latexpdf
mv build/latex/amrex.pdf source/
make clean
make SPHINXOPTS='-v -W --keep-going' PYTHON="python3" html
cd ../../

# copy sphinx to target location
cd build
cp -rp ../Docs/sphinx_documentation/build/html/* docs_html/
