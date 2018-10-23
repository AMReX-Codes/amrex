# AMREX_HOME defines the directory in which we will find all the BoxLib code
# If you set AMREX_HOME as an environment variable, this line will be ignored
AMREX_HOME ?= ../../../amrex

DEBUG     = FALSE
DEBUG     = TRUE
USE_MPI   = TRUE
USE_OMP   = FALSE
PROFILE   = FALSE
COMM_PROFILE  = FALSE
TRACE_PROFILE = FALSE
COMP      = gnu
#COMP      = Intel
DIM       = 3
PRECISION = DOUBLE
EBASE     = mfreadwrite

include ./Make.package
include $(AMREX_HOME)/Tools/GNUMake/Make.defs
include $(AMREX_HOME)/Src/Base/Make.package

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
