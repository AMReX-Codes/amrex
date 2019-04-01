# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line.
PYTHON       ?= python
SPHINXOPTS    =
SPHINXBUILD   = $(PYTHON) -msphinx
SPHINXPROJ    = amrex
SOURCEDIR     = source
BUILDDIR      = build

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help Makefile

clean:
	rm -rf build ../Doxygen/xml source/class source/file
	rm -f source/classlist.rst source/filelist.rst
	rm -rf source/*_files.rst

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)
