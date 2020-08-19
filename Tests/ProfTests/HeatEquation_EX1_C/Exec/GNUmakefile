# AMREX_HOME defines the directory in which we will find all the AMReX code.
AMREX_HOME = ../../../../

DEBUG         = FALSE
USE_MPI       = TRUE
USE_OMP       = FALSE
PROFILE       = TRUE
COMM_PROFILE  = TRUE
TRACE_PROFILE = TRUE
COMP          = gnu
DIM           = 2

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include ../Source/Make.package
VPATH_LOCATIONS  += ../Source

include $(AMREX_HOME)/Src/Base/Make.package

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
