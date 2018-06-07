# AMREX_HOME defines the directory in which we will find all the AMReX code.
# If you set AMREX_HOME as an environment variable, this line will be ignored
AMREX_HOME ?= ../../../../

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
