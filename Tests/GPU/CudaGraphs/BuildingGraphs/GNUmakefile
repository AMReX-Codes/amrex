AMREX_HOME ?= ../../../

DEBUG	= FALSE
#DEBUG	= TRUE

DIM	= 3

COMP    = gnu 

USE_MPI   = FALSE 
USE_OMP   = FALSE
USE_CUDA  = TRUE
TINY_PROFILE = TRUE

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include ./Make.package
include $(AMREX_HOME)/Src/Base/Make.package

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
