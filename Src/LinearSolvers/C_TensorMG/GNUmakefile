#
# $Id: GNUmakefile,v 1.8 1999-05-24 18:11:46 car Exp $
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

ifeq ($(KCC_VERSION),3.3)
CXXFLAGS += --one_instantiation_per_object
endif

all: $(optionsLib)

#
# Libraries to close against.
#
ifeq ($(COMP), KCC)
LibsToCloseAgainst := $(TOP)/lib/$(machineSuffix)/libmg$(DIM)d.a
LibsToCloseAgainst += $(TOP)/lib/$(machineSuffix)/libamr$(DIM)d.a
LibsToCloseAgainst += $(TOP)/lib/$(machineSuffix)/libbndry$(DIM)d.a
LibsToCloseAgainst += $(TOP)/lib/$(machineSuffix)/libbox$(DIM)d.a

$(optionsLib): $(LibsToCloseAgainst)

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

ifeq ($(MACHINE),T3E)
  #
  # Turn off warning about double precision constants.
  #
  FOPTF += -M 1110
  FDEBF += -M 1110
endif

include $(TOP)/mk/Make.rules
