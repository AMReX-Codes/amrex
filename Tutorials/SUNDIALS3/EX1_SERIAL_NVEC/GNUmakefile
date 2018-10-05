# AMREX_HOME defines the directory in which we will find all the BoxLib code
# If you set AMREX_HOME as an environment variable, this line will be ignored
AMREX_HOME ?= ../../../

DEBUG            = FALSE
USE_MPI          = TRUE #True should also work here
USE_OMP          = FALSE
USE_CUDA         = FALSE
COMP             = pgi # gnu should also work here
DIM              = 3
PRECISION        = DOUBLE
EBASE            = main

CVODE_LIB_DIR ?= ${CVODE_LIB}
USE_SUNDIALS3 = TRUE
USE_SUNDIALS3_INTERFACE = FALSE

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include ./Make.package
include $(AMREX_HOME)/Src/Base/Make.package

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
