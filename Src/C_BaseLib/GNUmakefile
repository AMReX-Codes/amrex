#
# $Id: GNUmakefile,v 1.10 1998-08-21 17:58:33 car Exp $
#

#
# Where is GUTHAMR?
#
GUTHAMR_HOME = ../IAMRAll

#
# PRECISION can be DOUBLE or FLOAT.
#
PRECISION = DOUBLE

#
# Debug can be TRUE or FALSE
#
DEBUG = TRUE

#
# DIM can be 1, 2, or 3.
#
DIM = 2
DIM = 3

#
# Your choice of C++ compiler: g++-2.5.8, OCC, CC, KCC, CFRONT
#
COMP = KCC

#
# The base name of the library we're building.
#
LBASE = box

#
# Should we use MPI?
#
USE_MPI = FALSE
USE_MPI = TRUE

PBOXLIB_HOME = ..
include ../mk/Make.defs Make.package

ifeq ($(USE_MPI),TRUE)
DEFINES += -DBL_USE_MPI
MPI_HOME = /usr/local/mpi
INCLUDE_LOCATIONS += $(MPI_HOME)/include
LIBRARY_LOCATIONS += $(MPI_HOME)/lib/alpha/ch_p4
LIBRARIES += -lmpi
endif

ifeq ($(KCC_VERSION),3.3)
CXXFLAGS+=--one_instantiation_per_object
endif


all: $(optionsLib)

#
# `install' does NOT try to build and install the HTML stuff.
# doc++ isn't widely distributed enough to assume that.
# For now you must invoke `install_html' separately.
#
install: all

include ../mk/Make.rules

CPPFLAGS += -DBL_USE_NEW_HFILES

#
# Temp stuff for doing html docs.
#
_docHeaders := $(sort $(wildcard *.H))

_htmlDir := html

_installDir = $(instHtmlDir)/boxlib

#
# This unconditionally installs stuff from ./html into _installDir.
#
install_html: html
	-if [ ! -d $(_installDir) ]; then \
              mkdir -p $(_installDir);    \
              chmod 755 $(_installDir);   \
        fi
	-cd $(_installDir); $(RM) *
	if [ -d $(_htmlDir) ]; then                     \
	    cd $(_htmlDir);                             \
	    for f in *; do                              \
	    	    $(INSTALL_DATA) $$f $(_installDir); \
	    done;                                       \
	fi

#
# This will make the HTML files in the directory html provided
# any of the relevant .H files have changed more recently than the
# directory itself.
#
html: $(_docHeaders)
	-if `doc++ -h > /dev/null` ; then                                 \
		$(RM) -r $(_htmlDir);                                     \
		mkdir -p $(_htmlDir);                                     \
		doc++ -f -j -B Banner.html -d $(_htmlDir) $(_docHeaders); \
	fi

BoxLib.tex: $(_docHeaders)
	-if `doc++ -h > /dev/null` ; then                                 \
	    doc++ -f -t -j -o BoxLib.tex -ep docxx_squish $(_docHeaders); \
	    sed '/[\]usepackage[{]docxx[}]/d' < BoxLib.tex > .junk;       \
	    mv .junk BoxLib.tex;                                          \
	fi

BoxLib.dvi: BoxLib.tex
	latex BoxLib.tex; latex BoxLib.tex

BoxLib.ps: BoxLib.dvi
	dvips BoxLib.dvi -o BoxLib.ps

#
# Install the PostScript file into GUTHMAR_HOME/ps
#
install_ps: BoxLib.ps
	-if [ ! -d $(GUTHAMR_HOME)/ps ]; then \
              mkdir -p $(GUTHAMR_HOME)/ps;    \
              chmod 755 $(GUTHAMR_HOME)/ps;   \
        fi
	$(INSTALL_DATA) BoxLib.ps $(GUTHAMR_HOME)/ps
