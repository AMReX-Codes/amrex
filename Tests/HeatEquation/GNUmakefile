# AMREX_HOME defines the directory in which we will find all the AMReX code.
AMREX_HOME ?= ../../../../amrex

USE_ASSERTION = TRUE

DEBUG        = FALSE
USE_MPI      = FALSE
USE_OMP      = FALSE
PRECISION	 = FLOAT
COMP         = gnu
DIM          = 3

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include Make.package
# VPATH_LOCATIONS  += ../Source
# INCLUDE_LOCATIONS += ../Source

include $(AMREX_HOME)/Src/Base/Make.package

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
