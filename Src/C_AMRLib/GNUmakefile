#
# $Id: GNUmakefile,v 1.10 1999-05-24 18:11:44 car Exp $
#
PBOXLIB_HOME = ..
TOP          = $(PBOXLIB_HOME)
#
# Variables for the user to set ...
#
PRECISION = DOUBLE
DEBUG     = TRUE
DIM       = 3
COMP      = KCC
USE_MPI   = TRUE
#
# The base name of the library we're building.
#
LBASE = amr
#
# Where libraries and include files will be installed.
#
INSTALL_ROOT = $(TOP)

include $(TOP)/mk/Make.defs
#
# Here's where we enable the accessibility of ArrayView.
#
USE_ARRAYVIEW = FALSE
#
# Use ArrayView on OSF1 ...
#
ifeq ($(MACHINE),OSF1)
USE_ARRAYVIEW = TRUE
endif

ifeq ($(USE_ARRAYVIEW),TRUE)
DEFINES += -DBL_USE_ARRAYVIEW
endif

include ./Make.package

INCLUDE_LOCATIONS += $(TOP)/include

ifeq ($(KCC_VERSION),3.3)
CXXFLAGS += --one_instantiation_per_object
endif

all: $(optionsLib)
#
# Libraries to close against.
#
ifeq ($(COMP), KCC)
LibsToCloseAgainst := $(TOP)/lib/$(machineSuffix)/libbndry$(DIM)d.a
LibsToCloseAgainst += $(TOP)/lib/$(machineSuffix)/libbox$(DIM)d.a

$(optionsLib): $(LibsToCloseAgainst)

endif
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

AmrLib.tex: $(_docHeaders)
	-if `doc++ -h > /dev/null` ; then                                 \
	    doc++ -f -t -j -o AmrLib.tex -ep docxx_squish $(_docHeaders); \
	    sed '/[\]usepackage[{]docxx[}]/d' < AmrLib.tex > .junk;       \
	    mv .junk AmrLib.tex;                                          \
	fi

AmrLib.dvi: AmrLib.tex
	latex AmrLib.tex; latex AmrLib.tex

AmrLib.ps: AmrLib.dvi
	dvips AmrLib.dvi -o AmrLib.ps
#
# Install the PostScript file into $(INSTALL_ROOT)/ps
#
install_ps: AmrLib.ps
	-if [ ! -d $(INSTALL_ROOT)/ps ]; then \
              mkdir -p $(INSTALL_ROOT)/ps;    \
              chmod 755 $(INSTALL_ROOT)/ps;   \
        fi
	$(INSTALL_DATA) AmrLib.ps $(INSTALL_ROOT)/ps
