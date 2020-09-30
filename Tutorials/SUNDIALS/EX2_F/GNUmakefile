# AMREX_HOME defines the directory in which we will find all the AMReX code
AMREX_HOME = ../../../

DEBUG     = FALSE
USE_MPI   = TRUE
USE_OMP   = FALSE
COMP      = gnu
DIM       = 3
PRECISION = DOUBLE
EBASE     = main

USE_SUNDIALS   = TRUE
USE_CVODE_LIBS = TRUE

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include ./Make.package
include $(AMREX_HOME)/Src/Base/Make.package

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
