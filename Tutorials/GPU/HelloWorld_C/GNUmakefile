AMREX_HOME ?= ../../../

DEBUG	= FALSE
#DEBUG	= TRUE

DIM	= 3

COMP    = pgi

USE_MPI   = TRUE 
USE_OMP   = FALSE
USE_CUDA  = TRUE

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include ./Make.package
include $(AMREX_HOME)/Src/Base/Make.package

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
