# AMREX_HOME defines the directory in which we will find all the AMReX code.
AMREX_HOME = ../../../../

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


DEBUG     = FALSE
USE_MPI   = FALSE
USE_OMP   = FALSE
COMP      = gnu
DIM       = 2

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include ../Source/Make.package
VPATH_LOCATIONS  += ../Source

include $(AMREX_HOME)/Src/Base/Make.package

ifeq ($(USE_CONDUIT),TRUE)
    include $(AMREX_HOME)/Src/Extern/Conduit/Make.package
endif

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
