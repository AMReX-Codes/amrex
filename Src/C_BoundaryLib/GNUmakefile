#
# $Id: GNUmakefile,v 1.2 1998-10-28 17:49:32 lijewski Exp $
#

PBOXLIB_HOME = ..

PRECISION = DOUBLE
DEBUG = TRUE
DIM = 3
COMP = KCC
USE_MPI = TRUE

#
# The base name of the library we're building.
#
LBASE = bndry

#
# Where libraries and include files will be installed.
#
INSTALL_ROOT = $(PBOXLIB_HOME)

include ../mk/Make.defs Make.package

INCLUDE_LOCATIONS += ../include

ifeq ($(USE_MPI), TRUE)
DEFINES += -DBL_USE_MPI
LIBRARIES += -lmpi
ifeq ($(MACHINE), OSF1)
MPI_HOME = /usr/local/mpi
INCLUDE_LOCATIONS += $(MPI_HOME)/include
LIBRARY_LOCATIONS += $(MPI_HOME)/lib/alpha/ch_p4
endif
endif

ifeq ($(KCC_VERSION),3.3)
CXXFLAGS += --one_instantiation_per_object
endif

#
# Libraries to close against.
#
ifeq ($(COMP), KCC)
LibsToCloseAgainst := ../lib/$(machineSuffix)/libbox$(DIM)d.a
endif

all: $(optionsLib)

#
# `install' does NOT try to build and install the HTML stuff.
# doc++ isn't widely distributed enough to assume that.
# For now you must invoke `install_html' separately.
#
install: all

include ../mk/Make.rules

#
# Temp stuff for doing html docs.
#
_docHeaders := $(sort $(wildcard *.H))

_htmlDir := html

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

BndryLib.tex: $(_docHeaders)
	-if `doc++ -h > /dev/null` ; then                                   \
	    doc++ -f -t -j -o BndryLib.tex -ep docxx_squish $(_docHeaders); \
	    sed '/[\]usepackage[{]docxx[}]/d' < BndryLib.tex > .junk;       \
	    mv .junk BndryLib.tex;                                          \
	fi

BndryLib.dvi: BndryLib.tex
	latex BndryLib.tex; latex BndryLib.tex

BndryLib.ps: BndryLib.dvi
	dvips BndryLib.dvi -o BndryLib.ps

#
# Install the PostScript file into $(INSTALL_ROOT)/ps
#
install_ps: BndryLib.ps
	-if [ ! -d $(INSTALL_ROOT)/ps ]; then \
              mkdir -p $(INSTALL_ROOT)/ps;    \
              chmod 755 $(INSTALL_ROOT)/ps;   \
        fi
	$(INSTALL_DATA) BndryLib.ps $(INSTALL_ROOT)/ps
