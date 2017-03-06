AMREX_HOME = ../../../amrex

DEBUG        = FALSE
USE_MPI      = TRUE
USE_OMP      = FALSE
USE_IPM      = FALSE
PROFILE      = TRUE
COMM_PROFILE = TRUE
USE_DUMPI    = FALSE
COMP         = g++
#COMP         = Intel
DIM          = 3
PRECISION    = DOUBLE
EBASE        = fbtest

CEXE_sources += MultiFabFillBoundaryTest.cpp

include $(AMREX_HOME)/Tools/GNUMake/Make.defs
include $(AMREX_HOME)/Src/Boundary/Make.package
include $(AMREX_HOME)/Src/Base/Make.package

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
