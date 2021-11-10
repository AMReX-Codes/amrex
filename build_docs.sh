#!/bin/bash
set -e # Exit with nonzero exit code if anything fails

# Doxygen
echo "Build the Doxygen documentation"
cd Docs/Doxygen
doxygen doxygen.conf &> doxygen.out
cd ../..

# sphinx
cd Docs/sphinx_documentation

echo "Build the Sphinx documentation for Amrex."
make PYTHON="python3" latexpdf
mv build/latex/amrex.pdf source/
make PYTHON="python3" html &> make_source_html.out
cd ../../

mkdir build
cd build
mkdir docs_html docs_xml

# add doxygen
mkdir -p docs_html/doxygen
cp -rp ../Docs/Doxygen/html/* docs_html/doxygen/
mkdir -p docs_xml/doxygen
cp -rp ../Docs/Doxygen/xml/* docs_xml/doxygen/
# add tagfile to allow other docs to interlink with amrex
cp ../Docs/Doxygen/amrex-doxygen-web.tag.xml docs_xml/doxygen/.

# add sphinx
cp -rp ../Docs/sphinx_documentation/build/html/* docs_html/

