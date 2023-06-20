# AMREX_HOME defines the directory in which we will find all the AMReX code.
# If you set AMREX_HOME as an environment variable, this line will be ignored
AMREX_HOME ?= ../../

DEBUG        = TRUE 
USE_MPI      = FALSE
MPI_THREAD_MULTIPLE = FALSE  # amrex.async_out=1 with more than 64 processes requires MPI_THREAD_MULTIPLE
USE_OMP      = FALSE
#USE_CUDA     = TRUE
COMP         = gnu
DIM          = 2

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include ./Make.package
include $(AMREX_HOME)/Src/Base/Make.package
include $(AMREX_HOME)/Src/AmrCore/Make.package
include $(AMREX_HOME)/Src/Boundary/Make.package

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
