AMREX_HOME ?= ../../..

DEBUG	= TRUE 

DIM	= 3

COMP    = hip

USE_CUDA  = FALSE 

USE_MPI   = FALSE
USE_OMP   = FALSE
USE_HIP   = TRUE

TINY_PROFILE = TRUE

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include ./Make.package
include $(AMREX_HOME)/Src/Base/Make.package

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
