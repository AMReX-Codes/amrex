AMREX_HOME ?= ../../../

DEBUG	= FALSE

DIM	= 3

COMP    = gcc

USE_CUDA  = TRUE
USE_ACC   = FALSE
USE_OMP_OFFLOAD = FALSE

USE_MPI   = FALSE
USE_OMP   = FALSE

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include ./Make.package
include $(AMREX_HOME)/Src/Base/Make.package

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
