#
# $Id: GNUmakefile,v 1.7 1999-05-24 18:11:45 car Exp $
#
PBOXLIB_HOME = ..

TOP = $(PBOXLIB_HOME)

PRECISION = DOUBLE
DEBUG = TRUE
DIM = 3
COMP = KCC
USE_MPI = TRUE
#
# The base name of the library we're building.
#
LBASE = mg
#
# Where libraries and include files will be installed.
#
INSTALL_ROOT = $(TOP)

include $(TOP)/mk/Make.defs Make.package

INCLUDE_LOCATIONS += $(TOP)/include

ifeq ($(KCC_VERSION),3.3)
CXXFLAGS += --one_instantiation_per_object
endif
#
# Libraries to close against.
#
ifeq ($(COMP), KCC)
LibsToCloseAgainst := $(TOP)/lib/$(machineSuffix)/libamr$(DIM)d.a
LibsToCloseAgainst += $(TOP)/lib/$(machineSuffix)/libbndry$(DIM)d.a
LibsToCloseAgainst += $(TOP)/lib/$(machineSuffix)/libbox$(DIM)d.a
endif

all: $(optionsLib)

include $(TOP)/mk/Make.rules

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

MgLib.tex: $(_docHeaders)
	-if `doc++ -h > /dev/null` ; then                                \
	    doc++ -f -t -j -o MgLib.tex -ep docxx_squish $(_docHeaders); \
	    sed '/[\]usepackage[{]docxx[}]/d' < MgLib.tex > .junk;       \
	    mv .junk MgLib.tex;                                          \
	fi

MgLib.dvi: MgLib.tex
	latex MgLib.tex; latex MgLib.tex

MgLib.ps: MgLib.dvi
	dvips MgLib.dvi -o MgLib.ps
#
# Install the PostScript file into $(INSTALL_ROOT)/ps
#
install_ps: MgLib.ps
	-if [ ! -d $(INSTALL_ROOT)/ps ]; then \
              mkdir -p $(INSTALL_ROOT)/ps;    \
              chmod 755 $(INSTALL_ROOT)/ps;   \
        fi
	$(INSTALL_DATA) MgLib.ps $(INSTALL_ROOT)/ps
