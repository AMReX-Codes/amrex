# AMREX_HOME defines the directory in which we will find all the BoxLib code
# If you set AMREX_HOME as an environment variable, this line will be ignored
AMREX_HOME ?= ../../../

DEBUG     = FALSE
USE_MPI   = TRUE
USE_OMP   = FALSE
USE_GPU   = TRUE
COMP      = gnu
DIM       = 3
PRECISION = DOUBLE
EBASE     = main

CVODE_LIB_DIR ?= ${CVODE_LIB}
SUNDIALS3_LIBS = -lsundials_nveccuda
USE_SUNDIALS3 = TRUE

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include ./Make.package
include $(AMREX_HOME)/Src/Base/Make.package

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
