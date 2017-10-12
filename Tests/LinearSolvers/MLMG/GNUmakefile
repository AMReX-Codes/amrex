AMREX_HOME ?= ../../../

DEBUG	= FALSE
DEBUG	= TRUE

DIM	= 3

COMP    = gnu

USE_MPI   = TRUE
USE_OMP   = TRUE

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include ./Make.package
include $(AMREX_HOME)/Src/Base/Make.package
include $(AMREX_HOME)/Src/Boundary/Make.package
include $(AMREX_HOME)/Src/LinearSolvers/MLMG/Make.package
include $(AMREX_HOME)/Src/LinearSolvers/C_CellMG/Make.package
include $(AMREX_HOME)/Src/AmrCore/Make.package

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
