AMREX_HOME ?= ../../../

DEBUG	= FALSE
#DEBUG	= TRUE

DIM	= 3

COMP    = gnu

USE_MPI   = TRUE
#USE_OMP   = TRUE
USE_OMP   = FALSE

USE_HYPRE = FALSE

TINY_PROFILE = TRUE

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include ./Make.package
include $(AMREX_HOME)/Src/Base/Make.package
include $(AMREX_HOME)/Src/Boundary/Make.package
include $(AMREX_HOME)/Src/LinearSolvers/MLMG/Make.package
include $(AMREX_HOME)/Src/LinearSolvers/C_CellMG/Make.package

ifeq ($(USE_HYPRE),TRUE)
    include $(AMREX_HOME)/Src/Extern/HYPRE/Make.package
endif

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
