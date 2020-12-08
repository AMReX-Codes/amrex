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
#
cd ../sphinx_tutorials
echo "Build the Sphinx documentation for the Amrex tutorials."
make PYTHON="python3" latexpdf &> make_tutorials_latex.out
mv build/latex/amrex.pdf source/
make PYTHON="python3" html &> make_tutorials_html.out
cd ../../

mkdir build
cd build
mkdir docs_html tutorials_html docs_xml

# add doxygen
mkdir -p docs_html/doxygen
cp -rp ../Docs/Doxygen/html/* docs_html/doxygen/
mkdir -p docs_xml/doxygen
cp -rp ../Docs/Doxygen/xml/* docs_xml/doxygen/

# add sphinx
cp -rp ../Docs/sphinx_documentation/build/html/* docs_html/
cp -rp ../Docs/sphinx_tutorials/build/html/* tutorials_html/
