# AMREX_HOME defines the directory in which we will find all the BoxLib code
# If you set AMREX_HOME as an environment variable, this line will be ignored
AMREX_HOME ?= ../..

DEBUG     = FALSE
USE_MPI   = FALSE
USE_OMP   = TRUE
PROFILE   = FALSE
COMP      = gnu
DIM       = 3
PRECISION = DOUBLE
EBASE     = Heat

include ./Make.package
include $(AMREX_HOME)/Tools/C_mk/Make.defs
include $(AMREX_HOME)/Src/C_BaseLib/Make.package

include $(AMREX_HOME)/Tools/C_mk/Make.rules
