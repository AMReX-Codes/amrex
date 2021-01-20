AMREX_HOME = ../../..

DEBUG	= TRUE
DEBUG	= FALSE

DIM	= 3

COMP    = gcc

USE_MPI   = FALSE
USE_OMP   = FALSE
USE_CUDA  = TRUE
USE_DPCPP = FALSE

## To test HIP indirect funciton, uncomment the following
# USE_HIP=TRUE
# HIP_INDIRECT_FUNCTION=TRUE

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include ./Make.package
include $(AMREX_HOME)/Src/Base/Make.package

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
