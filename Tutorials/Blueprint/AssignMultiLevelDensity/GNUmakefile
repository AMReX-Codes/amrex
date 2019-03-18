AMREX_HOME ?= ../../../

#
# This example requires installs of:
#  conduit (https://github.com/llnl/conduit)
#  ascent (https://github.com/Alpine-DAV/ascent)
#
# Example build:
#  make USE_CONDUIT=TRUE USE_ASCENT=TRUE \
#       CONDUIT_DIR={path/to/conduit/install} \
#       ASCENT_DIR={path/to/ascent/install}
#

DEBUG = FALSE

DIM	= 3

COMP    = gcc

TINY_PROFILE = FALSE
USE_PARTICLES = TRUE

PRECISION = DOUBLE

USE_MPI   = FALSE
USE_OMP   = FALSE

###################################################

EBASE     = main

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include ./Make.package
include $(AMREX_HOME)/Src/Base/Make.package
include $(AMREX_HOME)/Src/Boundary/Make.package
include $(AMREX_HOME)/Src/AmrCore/Make.package
include $(AMREX_HOME)/Src/Particle/Make.package

ifeq ($(USE_CONDUIT),TRUE)
    include $(AMREX_HOME)/Src/Extern/Conduit/Make.package
endif

include $(AMREX_HOME)/Tools/GNUMake/Make.rules

