#
# $Id: GNUmakefile,v 1.1 1997-09-15 19:40:05 lijewski Exp $
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

#
# TODO - remove this.
#
MAKEARCHFILE = ../IAMRAll/Make.OSF1

#
# Your choice of C++ compiler: g++-2.5.8, OCC, CC, KCC, CFRONT
#
COMP = KCC

#
# The base name of the library we're building.
#
LBASE = box

include $(GUTHAMR_HOME)/Make.defs Make.package

all: $(optionsLib)

#
# `install' does NOT try to build and install the HTML stuff.
# doc++ isn't widely distributed enough to assume that.
# For now you must invoke `install_html' separately.
#
install: all

include $(GUTHAMR_HOME)/Make.rules

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
