AMREX_HOME = ../../..

DEBUG	= TRUE
DEBUG	= FALSE

DIM	= 3

COMP    = gcc

USE_MPI   = FALSE
USE_OMP   = FALSE
USE_CUDA  = TRUE
USE_DPCPP = FALSE

## To test HIP indirect funciton, uncomment the following or run with `make -j8 USE_HIP=TRUE`
# USE_HIP=TRUE
ifdef USE_HIP
  HIP_INDIRECT_FUNCTION=$(USE_HIP)
endif

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include ./Make.package
include $(AMREX_HOME)/Src/Base/Make.package

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
