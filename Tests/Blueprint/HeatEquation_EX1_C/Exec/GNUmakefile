# AMREX_HOME defines the directory in which we will find all the AMReX code.
# If you set AMREX_HOME as an environment variable, this line will be ignored
AMREX_HOME ?= ../../../../

#
# This example requires installs of:
#  conduit (https://github.com/llnl/conduit)
#  ascent (https://github.com/Alpine-DAV/ascent)
#
# Example build:
#  make USE_CONDUIT=TRUE \
#       CONDUIT_DIR={path/to/conduit/install} \
#       ASCENT_DIR={path/to/ascent/install}
#


DEBUG     = FALSE
USE_MPI   = FALSE
USE_OMP   = FALSE
COMP      = gnu
DIM       = 2

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include ../Source/Make.package
VPATH_LOCATIONS  += ../Source

include $(AMREX_HOME)/Src/Base/Make.package

include $(AMREX_HOME)/Tools/GNUMake/Make.rules

ifeq ($(USE_CONDUIT), TRUE)
  # conduit + ascent
  INCLUDE_LOCATIONS += $(CONDUIT_DIR)/include
  INCLUDE_LOCATIONS += $(CONDUIT_DIR)/include/conduit
  VPATH_LOCATIONS   += $(CONDUIT_DIR)/include
  VPATH_LOCATIONS   += $(CONDUIT_DIR)/include/conduit
  LIBRARY_LOCATIONS += $(CONDUIT_DIR)/lib
  LIBRARIES += -lconduit -lconduit_blueprint -lconduit_relay
  
  INCLUDE_LOCATIONS += $(ASCENT_DIR)/include/ascent
  VPATH_LOCATIONS   += $(ASCENT_DIR)/include/ascent
  LIBRARY_LOCATIONS += $(ASCENT_DIR)/lib

ifeq ($(USE_MPI), TRUE)
  LIBRARIES += -lascent_mpi -lconduit_relay_mpi
else
  LIBRARIES += -lascent
endif

  DEFINES += -DUSE_CONDUIT
endif

