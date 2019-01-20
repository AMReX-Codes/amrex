# AMREX_HOME defines the directory in which we will find all the AMReX code.
# If you set AMREX_HOME as an environment variable, this line will be ignored
AMREX_HOME ?= ../../../../

DIM       = 3 

DEBUG     = FALSE
USE_MPI   = TRUE 
USE_OMP   = FALSE
USE_CUDA  = TRUE 

COMP      = pgi
USE_ACC   = TRUE

#COMP      = ibm
#USE_OMP_OFFLOAD = FALSE

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include ../Source_rewrite/Make.package
VPATH_LOCATIONS  += ../Source_rewrite

include $(AMREX_HOME)/Src/Base/Make.package

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
