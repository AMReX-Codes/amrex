AMREX_HOME = ../../../

DEBUG	= FALSE
DIM	= 3
COMP    = gcc

USE_MPI   = TRUE
USE_OMP   = FALSE
BL_NO_FORT = TRUE

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include ./Make.package
include $(AMREX_HOME)/Src/Base/Make.package
include $(AMREX_HOME)/Src/Boundary/Make.package
include $(AMREX_HOME)/Src/AmrCore/Make.package

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
