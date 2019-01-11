AMREX_HOME ?= ../../../

DEBUG	= FALSE

DIM	= 3

COMP    = ibm

USE_CUDA  = TRUE
USE_ACC   = FALSE

USE_MPI   = FALSE
USE_OMP   = FALSE
USE_OMP_OFFLOAD = TRUE

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include ./Make.package
include $(AMREX_HOME)/Src/Base/Make.package

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
