# AMREX_HOME defines the directory in which we will find all the BoxLib code
# If you set AMREX_HOME as an environment variable, this line will be ignored
AMREX_HOME ?= ../../../

DEBUG     = FALSE
USE_MPI   = FALSE
USE_OMP   = FALSE
COMP      = gnu
DIM       = 3
PRECISION = DOUBLE
EBASE     = ark_analytic

CVODE_LIB_DIR ?= ${CVODE_LIB}
USE_ARKODE_LIBS = TRUE
USE_SUNDIALS_3x4x = TRUE
SUNDIALS_VERSION = 4

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include ./Make.package
include $(AMREX_HOME)/Src/Base/Make.package

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
