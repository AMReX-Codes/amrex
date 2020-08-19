# AMREX_HOME defines the directory in which we will find all the AMReX code.
AMREX_HOME = ../../../..

DEBUG     = FALSE
USE_MPI   = TRUE 
USE_OMP   = FALSE
COMP      = intel 
DIM       = 2
THIRD_PARTY_PROFILING=CRAYPAT

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include ../Source/Make.package
VPATH_LOCATIONS  += ../Source

include $(AMREX_HOME)/Src/Base/Make.package

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
