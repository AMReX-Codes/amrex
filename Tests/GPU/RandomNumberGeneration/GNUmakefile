AMREX_HOME ?= ../../../

DEBUG	= FALSE

DIM	= 3

COMP    = hip

USE_CUDA  = FALSE 
USE_HIP   = TRUE
USE_ACC   = FALSE
USE_OMP_OFFLOAD = FALSE

USE_MPI   = FALSE
USE_OMP   = FALSE

TINY_PROFILE = TRUE

ifeq ($(USE_HIP),TRUE)
  INCLUDE_LOCATIONS += /opt/rocm/hiprand/include /opt/rocm/rocrand/include
  LIBRARY_LOCATIONS += /opt/rocm/hiprand/lib/    /opt/rocm/rocrand/lib/
  LIBRARIES += -lhiprand -lrocrand
endif

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include ./Make.package
include $(AMREX_HOME)/Src/Base/Make.package

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
