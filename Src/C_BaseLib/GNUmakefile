#
# $Id: GNUmakefile,v 1.21 2001-07-23 04:32:20 car Exp $
#

PBOXLIB_HOME = ..

TOP = $(PBOXLIB_HOME)

PRECISION = DOUBLE
DEBUG     = TRUE
DIM       = 3
COMP      = KCC
USE_MPI   = TRUE
#NAMESPACE = TRUE
NAMESPACE = FALSE

#
# The base name of the library we're building.
#
LBASE = box
#
# Where libraries and include files will be installed.
#
INSTALL_ROOT = $(PBOXLIB_HOME)

include $(TOP)/mk/Make.defs ./Make.package

ifeq ($(KCC_VERSION),3.3)
CXXFLAGS += --one_instantiation_per_object
endif

all: $(optionsLib)
#
# `install' does NOT try to build and install the HTML stuff.
# doc++ isn't widely distributed enough to assume that.
# For now you must invoke `install_html' separately.
#
install: all

include $(TOP)/mk/Make.rules

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
# Install the PostScript file into ./ps
#
install_ps: BoxLib.ps
	-if [ ! -d $(INSTALL_ROOT)/ps ]; then \
              mkdir -p $(INSTALL_ROOT)/ps;    \
              chmod 755 $(INSTALL_ROOT)/ps;   \
        fi
	$(INSTALL_DATA) BoxLib.ps $(INSTALL_ROOT)/ps

WRAPPERGEN=$(PBOXLIB_HOME)/util/wrappergen/wrappergen
BLMpi.cpp: BLMpi.w
	$(WRAPPERGEN) -c++ -o BLMpi.cpp -w BLMpi.w

clean::
	$(RM) -r html
