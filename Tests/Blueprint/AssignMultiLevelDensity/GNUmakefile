AMREX_HOME ?= ../../../

#
# This examples requires conduit (https://github.com/llnl/conduit)
# Example build:
#  make USE_CONDUIT=TRUE CONDUIT_DIR={path/to/conduit/install}
#

DEBUG = FALSE

DIM	= 3

COMP    = gcc

TINY_PROFILE = FALSE
USE_PARTICLES = TRUE

PRECISION = DOUBLE

USE_MPI   = TRUE
USE_OMP   = FALSE

###################################################

EBASE     = main

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include ./Make.package
include $(AMREX_HOME)/Src/Base/Make.package
include $(AMREX_HOME)/Src/Boundary/Make.package
include $(AMREX_HOME)/Src/AmrCore/Make.package
include $(AMREX_HOME)/Src/Particle/Make.package

include $(AMREX_HOME)/Tools/GNUMake/Make.rules

ifeq ($(USE_CONDUIT), TRUE)
  # conduit
  INCLUDE_LOCATIONS += $(CONDUIT_DIR)/include
  VPATH_LOCATIONS   += $(CONDUIT_DIR)/include
  LIBRARY_LOCATIONS += $(CONDUIT_DIR)/lib
  LIBRARIES += -lconduit -lconduit_blueprint -lconduit_relay

  DEFINES += -DUSE_CONDUIT
endif

