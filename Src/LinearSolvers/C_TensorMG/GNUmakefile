#
# $Id: GNUmakefile,v 1.2 1998-11-18 20:34:50 lijewski Exp $
#
PBOXLIB_HOME = ..

TOP = $(PBOXLIB_HOME)

PRECISION  = DOUBLE
DEBUG	   = FALSE
DEBUG	   = TRUE
PROFILE    = FALSE
DIM	   = 3
DIM	   = 2
COMP       = KCC
USE_MPI    = FALSE
USE_MPI    = TRUE
USE_NETCDF = FALSE
#
# The base name of the library we're building.
#
LBASE = mcmg
#
# Where libraries and include files will be installed.
#
INSTALL_ROOT = $(TOP)

include $(TOP)/mk/Make.defs ./Make.package

INCLUDE_LOCATIONS += $(TOP)/include

ifeq ($(USE_MPI), TRUE)
LIBRARIES += -lmpi
ifeq ($(MACHINE), OSF1)
INCLUDE_LOCATIONS += /usr/local/mpi/include
LIBRARY_LOCATIONS += /usr/local/mpi/lib/alpha/ch_p4
endif
endif

ifeq ($(KCC_VERSION),3.3)
CXXFLAGS += --one_instantiation_per_object
endif

#
# Libraries to close against.
#
ifeq ($(COMP), KCC)
LibsToCloseAgainst := $(TOP)/lib/$(machineSuffix)/libmg$(DIM)d.a
LibsToCloseAgainst += $(TOP)/lib/$(machineSuffix)/libamr$(DIM)d.a
LibsToCloseAgainst += $(TOP)/lib/$(machineSuffix)/libbndry$(DIM)d.a
LibsToCloseAgainst += $(TOP)/lib/$(machineSuffix)/libbox$(DIM)d.a
endif

ifeq ($(USE_NETCDF),TRUE)
LIBRARIES += /usr/people/stevens/bin/libnetcdf.a
INCLUDE_LOCATIONS += /usr/people/stevens/bin
endif
#
# For running 3rd only
#
3RD = 1
3RD =
ifdef 3RD
LDFLAGS += --link_command_prefix 3rd
LDFLAGS += -non_shared -v
endif

CXXFLAGS +=
CXXOPTF  +=
CXXDEBF  +=

CFLAGS +=
COPTF  +=
CDEBF  +=

FOPTF = -fast

all: $(optionsLib)

godzillaLink:
	rsh godzilla 'cd $(PWD); make $(executable)

include $(TOP)/mk/Make.rules
